import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Iterator, Callable, Any, Optional, List, Dict

import itertools
import pickle
import pandas as pd
import numpy as np
import scipy.special
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

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
    parser.add_argument('--mdsine_outdir', type=str, required=True)
    parser.add_argument('--mdsine_nomodule_outdir', type=str, required=True)
    parser.add_argument('--clv_elastic_outdir', type=str, required=True)
    parser.add_argument('--glv_elastic_outdir', type=str, required=True)
    parser.add_argument('--glv_ra_elastic_outdir', type=str, required=True)
    parser.add_argument('--glv_ra_ridge_outdir', type=str, required=True)
    parser.add_argument('--glv_ridge_outdir', type=str, required=True)
    parser.add_argument('--lra_elastic_outdir', type=str, required=True)
    parser.add_argument('--plot_dir', type=str, required=True)
    parser.add_argument('--sim_dt', type=float, required=False, default=0.01)
    parser.add_argument('--sim_max', type=float, required=False, default=1e20)
    parser.add_argument('--subsample_every', type=int, required=False, default=1)
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

    def evaluate_relative(self, pred: np.ndarray) -> pd.DataFrame:
        """Compute RMS error metric between prediction and truth (in relative abundance), one metric for each taxa."""
        truth = self.trajectory_subset(self.subject.times[0], self.subject.times[-1])
        rel_truth = truth / truth.sum(axis=0, keepdims=True)
        # rel_truth[rel_truth < lower_bound] = lower_bound

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
    mdsine2: Path
    mdsine2_nomodule: Path
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
        _require_file(self.mdsine2_nomodule)
        _require_file(self.clv_elastic)
        _require_file(self.glv_elastic)
        _require_file(self.glv_ra_elastic)
        _require_file(self.glv_ra_ridge)
        _require_file(self.glv_ridge)
        _require_file(self.lra_elastic)

    def mdsine2_fwsim(self, heldout: HoldoutData, sim_dt: float, sim_max: float, subsample_every: int = 1) -> np.ndarray:
        return forward_sim_mdsine2(data_path=self.mdsine2,
                                   recompute_cache=self.recompute_cache,
                                   heldout=heldout, sim_dt=sim_dt, sim_max=sim_max, init_limit_of_detection=1e5,
                                   subsample_every=subsample_every)

    def mdsine2_nomodule_fwsim(self, heldout: HoldoutData, sim_dt: float, sim_max: float, subsample_every: int = 1) -> np.ndarray:
        return forward_sim_mdsine2(data_path=self.mdsine2_nomodule,
                                   recompute_cache=self.recompute_cache,
                                   heldout=heldout, sim_dt=sim_dt, sim_max=sim_max, init_limit_of_detection=1e5,
                                   subsample_every=subsample_every)

    def clv_elastic_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray) -> np.ndarray:
        return forward_sim_clv(data_path=self.clv_elastic, recompute_cache=self.recompute_cache, x0=x0, u=u, t=t)

    def glv_elastic_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray, scale: float) -> np.ndarray:
        return forward_sim_glv(data_path=self.glv_elastic, recompute_cache=self.recompute_cache,
                               x0=x0, u=u, t=t, scale=scale, rel_abund=False, init_limit_of_detection=1e5)

    def glv_ra_elastic_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray, scale: float) -> np.ndarray:
        return forward_sim_glv(data_path=self.glv_ra_elastic, recompute_cache=self.recompute_cache,
                               x0=x0, u=u, t=t, scale=scale, rel_abund=True, init_limit_of_detection=1e5)

    def glv_ra_ridge_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray, scale: float) -> np.ndarray:
        return forward_sim_glv(data_path=self.glv_ra_ridge, recompute_cache=self.recompute_cache,
                               x0=x0, u=u, t=t, scale=scale, rel_abund=True, init_limit_of_detection=1e5)

    def glv_ridge_fwsim(self, x0: np.ndarray, u: np.ndarray, t: np.ndarray, scale: float) -> np.ndarray:
        return forward_sim_glv(data_path=self.glv_ridge, recompute_cache=self.recompute_cache,
                               x0=x0, u=u, t=t, scale=scale, rel_abund=False, init_limit_of_detection=1e5)

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
            mdsine2_nomodule=directories.mdsine2_nomodule / subject_id / "healthy" / "mcmc.pkl",
            clv_elastic=directories.clv_elastic / f'clv-{subject_idx}-model.pkl',
            glv_elastic=directories.glv_elastic / f'glv-elastic-net-{subject_idx}-model.pkl',
            glv_ra_elastic=directories.glv_ra_elastic / f'glv-ra-elastic-net-{subject_idx}-model.pkl',
            glv_ra_ridge=directories.glv_ra_ridge / f'glv-ra-ridge-{subject_idx}-model.pkl',
            glv_ridge=directories.glv_ridge / f'glv-ridge-{subject_idx}-model.pkl',
            lra_elastic=directories.lra_elastic / f'lra-{subject_idx}-model.pkl',
            recompute_cache=directories.recompute_cache,
        )


def evaluate_all(regression_inputs_dir: Path,
                 directories: HeldoutInferences,
                 complete_study: md2.Study,
                 sim_dt: float,
                 sim_max: float,
                 mdsine2_subsample_every: int = 1):
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
    for sidx, sid, inferences in retrieve_grouped_results(directories):
        x0, u, t = Y[sidx][0], U[sidx], T[sidx]

        heldout_data = HoldoutData(complete_study[sid], sidx)

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

        def add_relative_results(_method: str, _res: pd.DataFrame):
            for _, row in _res.iterrows():
                relative_df_entries.append({
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
            traj = np.median(
                inferences.mdsine2_fwsim(heldout_data, sim_dt, sim_max, subsample_every=mdsine2_subsample_every),
                axis=0
            )
            add_absolute_results(
                'MDSINE2',
                heldout_data.evaluate_absolute(traj)
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate MDSINE2 output: Holdout Subject {sid}.")

        try:
            add_absolute_results(
                'MDSINE2 (No Modules)',
                heldout_data.evaluate_absolute(
                    np.median(
                        inferences.mdsine2_nomodule_fwsim(heldout_data, sim_dt, sim_max, subsample_every=mdsine2_subsample_every),
                        axis=0
                    )
                )
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate MDSINE2 (nomodule) output: Holdout Subject {sid}.")

        try:
            add_absolute_results(
                'gLV (elastic net)',
                heldout_data.evaluate_absolute(inferences.glv_elastic_fwsim(x0, u, t, scale))
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate glv-elastic Net output: Holdout Subject {sid}.")

        try:
            add_absolute_results(
                'gLV (ridge)',
                heldout_data.evaluate_absolute(inferences.glv_ridge_fwsim(x0, u, t, scale))
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate glv-ridge output: Holdout Subject {sid}.")

        # Relative abundance
        try:
            add_relative_results(
                'MDSINE2',
                heldout_data.evaluate_relative(
                    np.median(
                        inferences.mdsine2_fwsim(heldout_data, sim_dt, sim_max, subsample_every=mdsine2_subsample_every),
                        axis=0
                    )
                )
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate MDSINE2 output (relabund): Holdout Subject {sid}.")

        try:
            add_relative_results(
                'MDSINE2 (No Modules)',
                heldout_data.evaluate_relative(
                    np.median(
                        inferences.mdsine2_nomodule_fwsim(heldout_data, sim_dt, sim_max, subsample_every=mdsine2_subsample_every),
                        axis=0
                    )
                )
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate MDSINE2 output (relabund): Holdout Subject {sid}.")

        try:
            add_relative_results(
                'cLV',
                heldout_data.evaluate_relative(inferences.clv_elastic_fwsim(x0, u, t))
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate clv output (relabund): Holdout Subject {sid}.")

        try:
            add_relative_results(
                'gLV-RA (elastic net)',
                heldout_data.evaluate_relative(inferences.glv_ra_elastic_fwsim(x0, u, t, scale))
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate glv-RA-elastic net output (relabund): Holdout Subject {sid}.")

        try:
            add_relative_results(
                'gLV-RA (ridge)',
                heldout_data.evaluate_relative(inferences.glv_ra_ridge_fwsim(x0, u, t, scale))
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate glv-RA-ridge output (relabund): Holdout Subject {sid}.")

        try:
            add_relative_results(
                'gLV (elastic net)',
                heldout_data.evaluate_relative(inferences.glv_elastic_fwsim(x0, u, t, scale))
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate glv-elastic net output (relabund): Holdout Subject {sid}.")

        try:
            add_relative_results(
                'gLV (ridge)',
                heldout_data.evaluate_relative(inferences.glv_ridge_fwsim(x0, u, t, scale))
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate glv-ridge output (relabund): Holdout Subject {sid}.")

        try:
            add_relative_results(
                'LRA',
                heldout_data.evaluate_relative(inferences.lra_elastic_fwsim(x0, u, t))
            )
        except FileNotFoundError:
            logger.error(f"Couldn't locate LRA output (relabund): Holdout Subject {sid}.")

    absolute_results = pd.DataFrame(absolute_df_entries)
    relative_results = pd.DataFrame(relative_df_entries)
    return absolute_results, relative_results


def make_boxplot(ax, df: pd.DataFrame,
                 method_order: List[str],
                 method_colors: Dict[str, np.ndarray],
                 lb: float,
                 ub: float,
                 xlabel: Optional[str] = None,
                 ylabel: Optional[str] = None):
    def agg_fn(_df):
        truth = _df['Truth'].to_numpy()
        pred = _df['Pred'].to_numpy()

        truth = np.where(truth < lb, lb, truth)
        truth = np.where(truth > ub, ub, truth)
        pred = np.where(pred < lb, lb, pred)
        pred = np.where(pred > ub, ub, pred)

        err = np.sqrt(np.mean(np.square(
            np.log10(pred) - np.log10(truth)
        )))  # Root mean squared?
        return pd.Series({'Error': err}, index=['Error'])

    df = df.groupby(['Method', 'HeldoutSubjectIdx', 'TaxonIdx']).apply(agg_fn).reset_index()
    df = df.loc[
        df['Method'].isin(method_order)
    ].sort_values(
        by='Method',
        key=lambda col: col.map({m: i for i, m in enumerate(method_order)})
    )

    sns.boxplot(
        data=df,
        ax=ax,
        x='Method',
        y='Error',
        showfliers=False,
        palette=method_colors,
        whis=(2.5, 97.5)
    )

    sns.stripplot(
        data=df,
        ax=ax,
        x='Method',
        y='Error',
        dodge=True,
        color='black',
        alpha=0.3,
        linewidth=1.0
    )

    def line_break(x: str):
        break_at_chars = ''.join(['(', ',', '/', ')'])
        return '\n'.join(re.split(f'[{break_at_chars}]', x))

    labels = [
        line_break(item.get_text())
        for item in ax.get_xticklabels()
    ]
    ax.set_xticklabels(labels)

    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=15)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=15)


def make_grouped_boxplot(abund_ax, error_ax,
                         df,
                         method_order: List[str],
                         method_colors: Dict[str, np.ndarray],
                         lb: float,
                         ub: float,
                         num_bins: int = 10,
                         error_ylabel: Optional[str] = None):
    def agg_fn(_df):
        truth = _df['Truth'].to_numpy()
        pred = _df['Pred'].to_numpy()

        truth = np.where(truth < lb, lb, truth)
        truth = np.where(truth > ub, ub, truth)
        pred = np.where(pred < lb, lb, pred)
        pred = np.where(pred > ub, ub, pred)

        err = np.sqrt(np.mean(np.square(
            np.log10(pred) - np.log10(truth)
        )))  # Root mean squared?
        return pd.Series({'Error': err}, index=['Error'])

    log_lb = np.log10(lb)
    df = df.assign(Bin=pd.qcut(
        np.log10(df['Truth']).replace([-np.inf], log_lb),
        q=num_bins,
        duplicates='drop'
    ))
    df_agg = df.groupby(['Method', 'HeldoutSubjectIdx', 'TaxonIdx', 'Bin']).apply(agg_fn).reset_index()

    #=====================================================================
    # Divide taxa based on abundance quantiles.
    # df = df.loc[
    #     df['Method'].isin(method_order)
    # ].assign(Bin=pd.qcut(df['TrueAbundMean'], q=num_bins))
    # log_lb = np.log10(lb)

    # ============ Render bin counts.
    def _aggregate_abundances_bin(_df):
        print(_df.head(1)['Bin'])
        bin = _df.head(1)['Bin'].item()
        count = sum(1 for _ in _df.groupby(['HeldoutSubjectIdx', 'TaxonIdx', 'TimePoint']))
        return pd.Series({
            'Left': bin.left if not np.isinf(bin.left) else log_lb,
            'Right': bin.right,
            'Count': count
        })
    bin_counts = df.groupby('Bin').apply(_aggregate_abundances_bin)
    widths = bin_counts['Right'] - bin_counts['Left']
    abund_ax.bar(
        x=bin_counts['Left'],
        height=1 / widths,
        align='edge',
        width=widths,
        linewidth=0.5,
        edgecolor='black',
        color='tab:blue'
    )
    abund_ax.set_xlabel('')
    abund_ax.set_ylabel('Bin Density')

    # ============= Render RMSE.
    def _bin_label(interval):
        left = interval.left
        right = interval.right
        if np.isinf(left):
            left = log_lb
        return '({:.1f},\n{:.1f}]'.format(left, right)
    df_agg['BinLabel'] = df_agg['Bin'].map(_bin_label)
    df_agg = df_agg.sort_values('Bin')
    supported_methods = set(pd.unique(df_agg['Method']))
    method_order = [m for m in method_order if m in supported_methods]
    sns.boxplot(
        data=df_agg,
        ax=error_ax,
        x='BinLabel',
        y='Error',
        hue='Method', hue_order=method_order, palette=method_colors,
        showfliers=False,
        whis=(2.5, 97.5)
    )
    error_ax.set_axisbelow(True)
    error_ax.yaxis.grid(True, 'major', linewidth=1, color='#e6e6e6')
    error_ax.set_xlabel('')
    if error_ylabel is not None:
        error_ax.set_ylabel(error_ylabel)

    abund_ax.legend([], [])
    error_ax.legend([], [])


def draw_method_legend(fig, method_order, method_colors, position: Tuple[float, float]):
    legend_elements = [
        Patch(facecolor=method_colors[method], edgecolor='black', linewidth=0.5, label=method)
        for method in method_order
    ]
    fig.legend(
        handles=legend_elements,
        loc='upper center', bbox_to_anchor=position,
        fancybox=True, ncol=len(method_order)
    )


def main():
    args = parse_args()
    plot_dir = Path(args.plot_dir)
    plot_dir.mkdir(exist_ok=True, parents=True)

    complete_study = md2.Study.load(args.study)
    directories = HeldoutInferences(
        mdsine2=Path(args.mdsine_outdir),
        mdsine2_nomodule=Path(args.mdsine_nomodule_outdir),
        clv_elastic=Path(args.clv_elastic_outdir),
        glv_elastic=Path(args.glv_elastic_outdir),
        glv_ra_elastic=Path(args.glv_ra_elastic_outdir),
        glv_ra_ridge=Path(args.glv_ra_ridge_outdir),
        glv_ridge=Path(args.glv_ridge_outdir),
        lra_elastic=Path(args.lra_elastic_outdir),
        recompute_cache=args.recompute_cache
    )

    # =================== Evaluate all errors.
    absolute_results, relative_results = evaluate_all(
        Path(args.regression_inputs_dir),
        directories,
        complete_study,
        args.sim_dt,
        args.sim_max,
        mdsine2_subsample_every=args.subsample_every
    )

    absolute_results.to_csv(plot_dir / "absolute_cv.tsv", sep='\t', index=False)
    relative_results.to_csv(plot_dir / "relative_cv.tsv", sep='\t', index=False)

    # ==================== Plot settings.
    methods = ['MDSINE2', 'MDSINE2 (No Modules)', 'cLV', 'LRA', 'gLV-RA (elastic net)', 'gLV-RA (ridge)', 'gLV (ridge)', 'gLV (elastic net)']
    palette_tab20 = sns.color_palette("tab10", len(methods))
    method_colors = {m: palette_tab20[i] for i, m in enumerate(methods)}

    # ==================== Evaluate and plot overall errors.
    fig, ax = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(14, 5),
        gridspec_kw={
            'width_ratios': [1, 2],
            'right': 0.92,
            'left': 0.08
        }
    )

    make_boxplot(ax[0], absolute_results, methods, method_colors, xlabel='Method', ylabel='RMSE (log Abs Abundance)', lb=1e-5, ub=1e40)
    make_boxplot(ax[1], relative_results, methods, method_colors, xlabel='Method', ylabel='RMSE (log Rel Abundance)', lb=1e-10, ub=1.0)
    fig.tight_layout()
    plt.savefig(plot_dir / "overall.pdf")
    plt.close(fig)

    # ==================== Evaluate and plot errors binned by log-abundances.
    fig, ax = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(20, 5),
        gridspec_kw={
            'height_ratios': [1, 2],
            'width_ratios': [1, 2],
            'wspace': 0.05,
            'right': 0.99,
            'left': 0.03,
            'bottom': 0.2,
            'top': 0.95
        }
    )

    make_grouped_boxplot(ax[0, 0], ax[1, 0], absolute_results, methods, method_colors, error_ylabel='RMSE (log Abs abundance)', lb=1e-5, ub=1e40)
    make_grouped_boxplot(ax[0, 1], ax[1, 1], relative_results, methods, method_colors, error_ylabel='RMSE (log Rel abundance)', lb=1e-10, ub=1.0)
    ax[0, 0].get_legend().remove()
    ax[1, 0].get_legend().remove()
    ax[0, 1].get_legend().remove()
    ax[1, 1].get_legend().remove()
    draw_method_legend(fig, methods, method_colors, position=(0.5, 0.1))
    fig.tight_layout()
    plt.savefig(plot_dir / "binned.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
