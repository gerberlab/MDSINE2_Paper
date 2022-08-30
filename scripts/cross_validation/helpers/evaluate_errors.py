import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Iterator, Callable, Any

import pickle
import pandas as pd
import numpy as np
import scipy.special
import matplotlib.pyplot as plt

import mdsine2 as md2
from mdsine2.names import STRNAMES
from mdsine2.logger import logger
from tqdm import tqdm

from lv_forward_sims import \
    add_limit_detection, \
    adjust_concentrations, \
    forward_sim_single_subj_clv, \
    forward_sim_single_subj_lra, \
    forward_sim_single_subj_glv


"""
MDSINE2 uses subject names, but clv code uses subject index.
"""
SUBJECT_IDS = ["2", "3", "4", "5"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=str, required=True)
    parser.add_argument('--regression_inputs_dir', type=str, required=True)
    parser.add_argument('--mdsine_outdir', type=str, required=True)
    parser.add_argument('--clv_elastic_outdir', type=str, required=True)
    parser.add_argument('--glv_elastic_outdir', type=str, required=True)
    parser.add_argument('--glv_ra_elastic_outdir', type=str, required=True)
    parser.add_argument('--glv_ra_ridge_outdir', type=str, required=True)
    parser.add_argument('--glv_ridge_outdir', type=str, required=True)
    parser.add_argument('--lra_elastic_outdir', type=str, required=True)
    parser.add_argument('--sim_dt', type=float, required=False, default=0.01)
    parser.add_argument('--sim_max', type=float, required=False, default=1e20)
    parser.add_argument('--recompute_cache', action='store_true')
    return parser.parse_args()


class HoldoutData:
    def __init__(self, subject: md2.Subject, subject_index: int):
        self.subject = subject
        self.subject_index = subject_index
        self.trajectories = self.subject.matrix()['abs']

    @property
    def initial_conditions(self) -> np.ndarray:
        return self.trajectories[:, 0]

    def trajectory_subset(self, start: float, end: float) -> np.ndarray:
        times = self.subject.times
        t_inside_range = (times >= start) & (times <= end)
        t_subset_indices, = np.where(t_inside_range)
        return self.trajectories[:, t_subset_indices]

    def evaluate_absolute(self, pred: np.ndarray, sim_max: float, lb: float = 1e5) -> float:
        """Compute RMS error metric between prediction and truth."""
        truth = self.trajectory_subset(self.subject.times[0], self.subject.times[-1])
        truth = np.where(truth < lb, lb, truth)
        pred = np.where(pred < lb, lb, pred)
        pred = np.where(pred > sim_max, sim_max, pred)
        if pred.shape != truth.shape:
            raise ValueError(f"truth shape ({truth.shape}) does not match pred shape ({pred.shape})")
        return np.sqrt(np.mean(np.square(pred - truth)))  # RMS

    def evaluate_relative(self, pred: np.ndarray) -> float:
        raise NotImplementedError()


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
def forward_sim_mdsine2(data_path: Path, heldout: HoldoutData, sim_dt: float, sim_max: float) -> np.ndarray:
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

            perts.append(mcmc.graph.perturbations[subj_pert.name].get_trace_from_disk())
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
        start_day=heldout.subject.times[0],
        sim_max=sim_max,
        perturbation_starts=pert_starts,
        perturbation_ends=pert_ends
    )

    n_samples = growth.shape[0]
    n_taxa = growth.shape[1]
    pred_matrix = np.zeros(shape=(n_samples, n_taxa, len(times)))
    for sample_idx in tqdm(range(n_samples), desc="MDSINE2 fwsim"):
        dyn.growth = growth[sample_idx]
        dyn.interactions = interactions[sample_idx]
        if perts is not None:
            dyn.perturbations = [pert[sample_idx] for pert in perts]

        init = heldout.initial_conditions
        if len(init.shape) == 1:
            init = init.reshape(-1, 1)
        x = md2.integrate(dynamics=dyn, initial_conditions=init,
                          dt=sim_dt, n_days=times[-1] + sim_dt, subsample=True, times=times)
        pred_matrix[sample_idx] = x['X']
    return pred_matrix


@cached_forward_simulation
def forward_sim_clv(data_path: Path,
                    x0: np.ndarray,
                    u: np.ndarray,
                    t: np.ndarray,
                    pseudo_count: int=0) -> np.ndarray:
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
                    rel_abund: bool) -> np.ndarray:
    logger.info(f"Evaluating gLV simulation using output ({data_path})")
    with open(data_path, "rb") as f:
        model = pickle.load(f)
        A, g, B = model.get_params()

    x0 = np.log(x0)
    if rel_abund:
        # Normalize (in log-scale)
        x0 = x0 - scipy.special.logsumexp(x0)
    else:
        """
        gLV inference is run with scaling (for numerical precision). 
        This is the inverse transformation!
        """
        A = A * scale
    return forward_sim_single_subj_glv(A, g, B, x0, u, t, rel_abund=rel_abund).transpose(1, 0)


@dataclass
class HeldoutInferences:
    mdsine2: Path
    clv_elastic: Path
    glv_elastic: Path
    glv_ra_elastic: Path
    glv_ra_ridge: Path
    glv_ridge: Path
    lra_elastic: Path
    recompute_cache: bool

    def __post_init__(self):
        def _require_file(path: Path):
            if not path.exists:
                raise FileNotFoundError(f"{path} not found!")

        _require_file(self.mdsine2)
        _require_file(self.clv_elastic)
        _require_file(self.glv_elastic)
        _require_file(self.glv_ra_elastic)
        _require_file(self.glv_ra_ridge)
        _require_file(self.glv_ridge)
        _require_file(self.lra_elastic)

    def mdsine2_fwsim(self, heldout: HoldoutData, sim_dt: float, sim_max: float) -> np.ndarray:
        return forward_sim_mdsine2(data_path=self.mdsine2, recompute_cache=self.recompute_cache, heldout=heldout, sim_dt=sim_dt, sim_max=sim_max)

    def clv_elastic_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray) -> np.ndarray:
        return forward_sim_clv(data_path=self.clv_elastic, recompute_cache=self.recompute_cache, x0=x0, u=u, t=t)

    def glv_elastic_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray, scale: float) -> np.ndarray:
        return forward_sim_glv(data_path=self.glv_elastic, recompute_cache=self.recompute_cache, x0=x0, u=u, t=t, scale=scale, rel_abund=False)

    def glv_ra_elastic_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray, scale: float) -> np.ndarray:
        return forward_sim_glv(data_path=self.glv_ra_elastic, recompute_cache=self.recompute_cache, x0=x0, u=u, t=t, scale=scale, rel_abund=True)

    def glv_ra_ridge_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray, scale: float) -> np.ndarray:
        return forward_sim_glv(data_path=self.glv_ra_ridge, recompute_cache=self.recompute_cache, x0=x0, u=u, t=t, scale=scale, rel_abund=True)

    def glv_ridge_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray, scale: float) -> np.ndarray:
        return forward_sim_glv(data_path=self.glv_ridge, recompute_cache=self.recompute_cache, x0=x0, u=u, t=t, scale=scale, rel_abund=False)

    def lra_elastic_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray) -> np.ndarray:
        return forward_sim_lra(data_path=self.lra_elastic, recompute_cache=self.recompute_cache, x0=x0, u=u, t=t)


def retrieve_grouped_results(directories: HeldoutInferences) -> Iterator[Tuple[int, str, HeldoutInferences]]:
    """
    Group the mdsine2/clv/glv etc result paths by heldout subject.

    :return: A (subject_id) -> (mdsine2, clv-elastic, glv-elastic, glv-ra-elastic, glv-ra-ridge, glv-ridge
    """
    for subject_idx, subject_id in enumerate(SUBJECT_IDS):
        yield subject_idx, subject_id, HeldoutInferences(
            mdsine2=directories.mdsine2 / subject_id / "healthy" / "mcmc.pkl",
            clv_elastic=directories.clv_elastic / f'clv-{subject_idx}-model.pkl',
            glv_elastic=directories.glv_elastic / f'glv-elastic-net-{subject_idx}-model.pkl',
            glv_ra_elastic=directories.glv_ra_elastic / f'glv-ra-elastic-net-{subject_idx}-model.pkl',
            glv_ra_ridge=directories.glv_ra_ridge / f'glv-ra-ridge-{subject_idx}-model.pkl',
            glv_ridge=directories.glv_ridge / f'glv-ridge-{subject_idx}-model.pkl',
            lra_elastic=directories.lra_elastic / f'lra-{subject_idx}-model.pkl',
            recompute_cache=directories.recompute_cache
        )


def evaluate_all(regression_inputs_dir: Path,
                 directories: HeldoutInferences,
                 complete_study: md2.Study,
                 sim_dt: float,
                 sim_max: float):
    # =========== Load regression inputs
    with open(regression_inputs_dir / "Y.pkl", "rb") as f:
        Y = pickle.load(f)
    with open(regression_inputs_dir / "U.pkl", "rb") as f:
        U = pickle.load(f)
    with open(regression_inputs_dir / "T.pkl", "rb") as f:
        T = pickle.load(f)

    # Include the limit of detection value
    Y = add_limit_detection(Y, 1e5)
    Y_adj, scale = adjust_concentrations(Y)
    logger.debug(f"Regression input rescaling value: {scale}")

    # =========== Load evaluations.
    absolute_df_entries = []
    relative_df_entries = []
    for sidx, sid, inferences in retrieve_grouped_results(directories):
        x0, u, t = Y[sidx][0], U[sidx], T[sidx]
        heldout_data = HoldoutData(complete_study[sid], sidx)

        def add_absolute_entry(_method, _err):
            absolute_df_entries.append({
                'HeldoutSubjectId': sid,
                'HeldoutSubjectIdx': sidx,
                'Method': _method,
                'Error': _err
            })

        def add_relative_entry(_method, _err):
            relative_df_entries.append({
                'HeldoutSubjectId': sid,
                'HeldoutSubjectIdx': sidx,
                'Method': _method,
                'Error': _err
            })

        # Absolute abundance
        add_absolute_entry(
            'MDSINE2',
            heldout_data.evaluate_absolute(np.median(inferences.mdsine2_fwsim(heldout_data, sim_dt, sim_max), axis=0), sim_max)
        )
        add_absolute_entry(
            'gLV-elastic net',
            heldout_data.evaluate_absolute(inferences.glv_elastic_fwsim(x0, u, t, scale), sim_max)
        )
        add_absolute_entry(
            'gLV-ridge',
            heldout_data.evaluate_absolute(inferences.glv_ridge_fwsim(x0, u, t, scale), sim_max)
        )

        # Relative abundance
        add_relative_entry(
            'MDSINE2',
            heldout_data.evaluate_relative(np.median(inferences.mdsine2_fwsim(heldout_data, sim_dt, sim_max), axis=0))
        )
        add_relative_entry(
            'cLV',
            heldout_data.evaluate_relative(inferences.clv_elastic_fwsim(x0, u, t))
        )
        add_relative_entry(
            'gLV-elastic net',
            heldout_data.evaluate_relative(inferences.glv_ra_elastic_fwsim(x0, u, t, scale))
        )
        add_relative_entry(
            'gLV-ridge',
            heldout_data.evaluate_relative(inferences.glv_ra_ridge_fwsim(x0, u, t, scale))
        )
        add_relative_entry(
            'LRA',
            heldout_data.evaluate_relative(inferences.lra_elastic_fwsim(x0, u, t))
        )

    absolute_results = pd.DataFrame(absolute_df_entries)
    relative_results = pd.DataFrame(relative_df_entries)
    return absolute_results, relative_results


def main():
    args = parse_args()

    complete_study = md2.Study.load(args.study)
    directories = HeldoutInferences(
        mdsine2=Path(args.mdsine_outdir),
        clv_elastic=Path(args.clv_elastic_outdir),
        glv_elastic=Path(args.glv_elastic_outdir),
        glv_ra_elastic=Path(args.glv_ra_elastic_outdir),
        glv_ra_ridge=Path(args.glv_ra_ridge_outdir),
        glv_ridge=Path(args.glv_ridge_outdir),
        lra_elastic=Path(args.lra_elastic_outdir),
        recompute_cache=args.recompute_cache
    )

    absolute_results, relative_results = evaluate_all(
        Path(args.regression_inputs_dir),
        directories,
        complete_study,
        args.sim_dt,
        args.sim_max
    )

    print(absolute_results)
    print(relative_results)

    # fig, ax = plt.subplots(
    #     nrows=1,
    #     ncols=4,
    #     figsize=(22, 4.5),
    #     gridspec_kw={
    #         'width_ratios': [1, 1, 2, 2]
    #     }
    # )


if __name__ == "__main__":
    main()