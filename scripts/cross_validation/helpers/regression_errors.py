import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Iterator, Callable, Any

import itertools
import pickle
import pandas as pd
import numpy as np
import scipy.special
import matplotlib

import mdsine2 as md2
from mdsine2.names import STRNAMES
from mdsine2.logger import logger
from tqdm import tqdm

from lv_forward_sims import \
    adjust_concentrations, \
    forward_sim_single_subj_clv, \
    forward_sim_single_subj_lra, \
    forward_sim_single_subj_glv

# Make font editable in AI
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

"""
MDSINE2 uses subject names, but clv code uses subject index.
"""
SUBJECT_IDS = ["2", "3", "4", "5"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=str, required=True)
    parser.add_argument('--regression_inputs_dir', type=str, required=True)
    parser.add_argument('--regression_outdir', type=str, required=True)
    parser.add_argument('--out_dir', type=str, required=True)
    parser.add_argument('--model_name', type=str, required=True)
    parser.add_argument('--recompute_cache', action='store_true')
    return parser.parse_args()


class HoldoutData:
    def __init__(self, subject: md2.Subject, subject_index: int):
        self.subject = subject
        self.subject_index = subject_index
        self.trajectories = self.subject.matrix()['abs']

    @property
    def normalized_trajectories(self) -> np.ndarray:
        return self.trajectories / self.trajectories.sum(axis=0, keepdims=True)

    def initial_conditions(self, lower_bound: float) -> np.ndarray:
        x2 = np.copy(self.trajectories[:, 0])
        x2[x2 < lower_bound] = lower_bound
        return x2

    def trajectory_subset(self, start: float, end: float) -> np.ndarray:
        """
        A (n_taxa) x (n_timepoints) matrix, confined to the time interval [start, end] (inclusive).
        """
        times = self.subject.times
        t_inside_range = (times >= start) & (times <= end)
        t_subset_indices, = np.where(t_inside_range)
        trajs = np.copy(self.trajectories[:, t_subset_indices])
        return trajs

    def evaluate_absolute(self, pred: np.ndarray) -> pd.DataFrame:
        """Compute RMS error metric between prediction and truth, one metric for each taxa."""
        truth = self.trajectory_subset(self.subject.times[0], self.subject.times[-1])
        # truth = np.where(truth < lower_bound, lower_bound, truth)
        # truth = np.where(truth > upper_bound, upper_bound, truth)
        #
        # pred = np.where(pred < lower_bound, lower_bound, pred)
        # pred = np.where(pred > upper_bound, upper_bound, pred)
        if pred.shape != truth.shape:
            raise ValueError(f"truth shape ({truth.shape}) does not match pred shape ({pred.shape})")

        # truth = np.log10(truth)
        # pred = np.log10(pred)
        entries = []
        for taxa_idx, time_idx in itertools.product(range(pred.shape[0]), range(pred.shape[1])):
            entries.append({
                'taxa': taxa_idx,
                'time': time_idx,
                'truth': truth[taxa_idx, time_idx],
                'pred': pred[taxa_idx, time_idx]
            })
        return pd.DataFrame(entries)

    def evaluate_relative(self, pred: np.ndarray, eps: float = 1e-10) -> pd.DataFrame:
        """Compute RMS error metric between prediction and truth (in relative abundance), one metric for each taxa."""
        truth = self.trajectory_subset(self.subject.times[0], self.subject.times[-1])
        rel_truth = truth / truth.sum(axis=0, keepdims=True)
        # rel_truth[rel_truth < lower_bound] = lower_bound

        pred = pred + eps
        rel_pred = pred / pred.sum(axis=0, keepdims=True)
        if rel_pred.shape != rel_truth.shape:
            raise ValueError(f"truth shape ({rel_truth.shape}) does not match pred shape ({rel_pred.shape})")

        # rel_pred[rel_pred < lower_bound] = lower_bound
        # rel_truth = np.log10(rel_truth)
        # rel_pred = np.log10(rel_pred)
        entries = []
        for taxa_idx, time_idx in itertools.product(range(rel_pred.shape[0]), range(rel_pred.shape[1])):
            entries.append({
                'taxa': taxa_idx,
                'time': time_idx,
                'truth': rel_truth[taxa_idx, time_idx],
                'pred': rel_pred[taxa_idx, time_idx]
            })
        return pd.DataFrame(entries)


def cached_forward_simulation(fwsim_fn: Callable[[Any], np.ndarray]):
    def _wrapper(*args, **kwargs):
        if 'data_path' not in kwargs:
            raise RuntimeError(f"function `{fwsim_fn.__name__}` should be called with a `data_path` kwarg.")

        if 'recompute_cache' in kwargs:
            recompute = 'recompute_cache' in kwargs and kwargs['recompute_cache'] == True
            del kwargs['recompute_cache']
        else:
            recompute = False

        data_path = kwargs['data_path']
        fwsim_path = data_path.with_suffix('.fwsim.npy')
        if recompute or (not fwsim_path.exists()):
            fwsim = fwsim_fn(*args, **kwargs)
            np.save(str(fwsim_path), fwsim)
            logger.debug(f"Saved forward simulation to {fwsim_path}")
            return fwsim
        else:
            logger.debug(f"Loaded previously evaluated forward simulation {fwsim_path}")
            return np.load(str(fwsim_path))
    return _wrapper


@cached_forward_simulation
def forward_sim_mdsine2(data_path: Path, heldout: HoldoutData, sim_dt: float, sim_max: float, init_limit_of_detection: float, subsample_every: int) -> np.ndarray:
    logger.info(f"Evaluating forward simulation using mdsine2 MCMC samples ({data_path})")
    mcmc = md2.BaseMCMC.load(str(data_path))

    growth = mcmc.graph[STRNAMES.GROWTH_VALUE].get_trace_from_disk()
    self_interactions = mcmc.graph[STRNAMES.SELF_INTERACTION_VALUE].get_trace_from_disk()
    interactions = mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk()
    interactions[np.isnan(interactions)] = 0
    self_interactions = -np.absolute(self_interactions)
    for i in range(self_interactions.shape[1]):
        interactions[:, i, i] = self_interactions[:, i]

    perts = []
    pert_starts = []
    pert_ends = []
    if heldout.subject.perturbations is not None:
        assert mcmc.graph.perturbations is not None
        for subj_pert in heldout.subject.perturbations:
            if subj_pert.name not in mcmc.graph.perturbations:
                raise KeyError(f"Heldout subject ({heldout.subject.name}) has perturbation `{subj_pert.name}`, "
                               f"but learned model does not.")

            pert_values = mcmc.graph.perturbations[subj_pert.name].get_trace_from_disk()
            pert_values[np.isnan(pert_values)] = 0.
            perts.append(pert_values)
            pert_starts.append(subj_pert.starts[heldout.subject.name])
            pert_ends.append(subj_pert.ends[heldout.subject.name])

        """
        Sort the perturbations in start order
        """
        sorted_order = np.argsort(pert_starts)
        perts = [perts[i] for i in sorted_order]
        pert_starts = [pert_starts[i] for i in sorted_order]
        pert_ends = [pert_ends[i] for i in sorted_order]
    else:
        perts = None
        pert_starts = None
        pert_ends = None

    times = heldout.subject.times
    dyn = md2.model.gLVDynamicsSingleClustering(
        growth=None,
        interactions=None,
        start_day=times[0],
        sim_max=sim_max,
        perturbation_starts=pert_starts,
        perturbation_ends=pert_ends
    )

    n_samples = growth.shape[0]
    n_taxa = growth.shape[1]

    gibbs_indices = list(range(0, n_samples, subsample_every))
    pred_matrix = np.empty(shape=(len(gibbs_indices), n_taxa, len(times)))

    for pred_idx, sample_idx in tqdm(enumerate(gibbs_indices), desc="MDSINE2 fwsim", total=len(pred_matrix)):
        dyn.growth = growth[sample_idx]
        dyn.interactions = interactions[sample_idx]
        if perts is not None:
            dyn.perturbations = [pert[sample_idx] for pert in perts]

        init = heldout.initial_conditions(lower_bound=init_limit_of_detection)
        if len(init.shape) == 1:
            init = init.reshape(-1, 1)
        x = md2.integrate(dynamics=dyn, initial_conditions=init,
                          dt=sim_dt, n_days=times[-1] + sim_dt, subsample=True, times=times)
        pred_matrix[pred_idx] = x['X']
    return pred_matrix


@cached_forward_simulation
def forward_sim_clv(data_path: Path,
                    x0: np.ndarray,
                    u: np.ndarray,
                    t: np.ndarray,
                    pseudo_count: int=1e-6) -> np.ndarray:
    logger.info(f"Evaluating cLV simulation using output ({data_path})")
    with open(data_path, "rb") as f:
        model = pickle.load(f)
        A, g, B = model.get_params()
        denom = model.denom

    x0 = x0 / np.sum(x0)  # normalize.
    return forward_sim_single_subj_clv(A, g, B, x0, u, t, denom, pc=pseudo_count).transpose(1, 0)


@cached_forward_simulation
def forward_sim_lra(data_path: Path,
                    x0: np.ndarray,
                    u: np.ndarray,
                    t: np.ndarray) -> np.ndarray:
    logger.info(f"Evaluating LRA simulation using output ({data_path})")
    with open(data_path, "rb") as f:
        model = pickle.load(f)
        A, g, B = model.get_params()

    x0 = x0 / np.sum(x0)  # normalize.
    return forward_sim_single_subj_lra(A, g, B, x0, u, t).transpose(1, 0)


@cached_forward_simulation
def forward_sim_glv(data_path: Path,
                    x0: np.ndarray,
                    u: np.ndarray,
                    t: np.ndarray,
                    scale: float,
                    init_limit_of_detection: float,
                    rel_abund: bool) -> np.ndarray:
    logger.info(f"Evaluating gLV simulation using output ({data_path})")
    with open(data_path, "rb") as f:
        model = pickle.load(f)
        A, g, B = model.get_params()

    x0 = np.copy(x0)
    x0[x0 < init_limit_of_detection] = init_limit_of_detection
    if rel_abund:
        # Normalize (in log-scale)
        x0 = np.log(x0)
        x0 = x0 - scipy.special.logsumexp(x0)
    else:
        """
        gLV inference is run with scaling (for numerical precision). 
        This is the inverse transformation!
        """
        # A = A * scale
        x0 = np.log(x0)

    # Include the limit of detection value
    return forward_sim_single_subj_glv(A, g, B, x0, u, t, rel_abund=rel_abund).transpose(1, 0)


@dataclass
class HeldoutInferences:
    regression_output: Path
    recompute_cache: bool

    def __post_init__(self):
        def _require_file(path: Path):
            if not path.exists:
                raise FileNotFoundError(f"{path} not found!")

        _require_file(self.regression_output)

    def glv_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray) -> np.ndarray:
        return forward_sim_glv(data_path=self.regression_output, recompute_cache=self.recompute_cache, x0=x0, u=u, t=t)


def evaluate_all(regression_inputs_dir: Path,
                 regression_outdir: Path,
                 model_name: str,
                 complete_study: md2.Study,
                 recompute_cache: bool):
    # =========== Load regression inputs
    with open(regression_inputs_dir / "Y.pkl", "rb") as f:
        Y = pickle.load(f)
    with open(regression_inputs_dir / "U.pkl", "rb") as f:
        U = pickle.load(f)
    with open(regression_inputs_dir / "T.pkl", "rb") as f:
        T = pickle.load(f)

    _, scale = adjust_concentrations(Y)
    logger.debug(f"Regression input rescaling value: {scale}")

    # =========== Load evaluations.
    absolute_df_entries = []
    relative_df_entries = []
    for sidx, sid in enumerate(SUBJECT_IDS):
        pkl_path = regression_outdir / f"{model_name}-{sidx}-model.pkl"
        heldout_data = HoldoutData(complete_study[sid], sidx)

        x0, u, t = Y[sidx][0], U[sidx], T[sidx]
        fwsim = forward_sim_glv(data_path=pkl_path, recompute_cache=recompute_cache,
                                x0=x0, u=u, t=t, scale=scale, rel_abund=False, init_limit_of_detection=1e5)

        def add_absolute_results(_method: str, _res: pd.DataFrame):
            for _, row in _res.iterrows():
                absolute_df_entries.append({
                    'HeldoutSubjectId': sid,
                    'HeldoutSubjectIdx': sidx,
                    'Method': _method,
                    'TaxonIdx': row['taxa'],
                    'TimePoint': row['time'],
                    'Truth': row['truth'],
                    'Pred': row['pred']
                })

        # Absolute abundance
        try:
            add_absolute_results(
                'gLV (elastic net)',
                heldout_data.evaluate_absolute(fwsim)
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate glv-elastic Net output: Holdout Subject {sid}.")

    absolute_results = pd.DataFrame(absolute_df_entries)
    relative_results = pd.DataFrame(relative_df_entries)
    return absolute_results, relative_results


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    complete_study = md2.Study.load(args.study)

    # =================== Evaluate all errors.
    absolute_results, relative_results = evaluate_all(
        Path(args.regression_inputs_dir),
        Path(args.regression_outdir),
        model_name=args.model_name,
        complete_study=complete_study,
        recompute_cache=args.recompute_cache
    )

    absolute_results.to_csv(out_dir / "absolute_cv.tsv", sep='\t', index=False)
    relative_results.to_csv(out_dir / "relative_cv.tsv", sep='\t', index=False)


if __name__ == "__main__":
    main()
