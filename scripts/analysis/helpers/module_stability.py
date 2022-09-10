import argparse
import itertools
from pathlib import Path
from typing import Tuple, Iterator

import numpy as np
import pandas as pd

import mdsine2 as md2
from mdsine2 import Clustering
from mdsine2.base import _Cluster
from mdsine2.names import STRNAMES


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--fixed-cluster-mcmc-path', '-f', type=str, dest='mcmc_path',
                        required=True,
                        help='Path of saved MDSINE2.BaseMCMC chain (fixed-clustering inference)')
    parser.add_argument('--study', '-s', dest='study', type=str, required=True,
                        help="The path to the relevant Study object containing the input data (subjects, taxa).")
    parser.add_argument('--module-remove-idx', '-m', dest='module_remove_idx', type=int,
                        help='Specify which module to remove, specified by index (zero-indexed)')

    parser.add_argument('--out-path', '-o', dest='out_path', required=True, type=str)

    parser.add_argument('--seed', required=False, type=int, default=31415)
    parser.add_argument('--simulation-dt', '-dt', type=float, dest='dt', required=False,
                        help='Timesteps we go in during forward simulation', default=0.01)
    parser.add_argument('--sim-max', dest='sim_max', type=float, required=False,
                        help='Maximum value of abundance.', default=1e20)
    return parser.parse_args()


def main():
    args = parse_args()

    np.random.seed(args.seed)

    study = md2.Study.load(args.study)
    mcmc = md2.BaseMCMC.load(args.mcmc_path)
    module_idx_to_remove = args.module_remove_idx

    modules: Clustering = mcmc.graph[STRNAMES.CLUSTERING_OBJ]
    module_to_remove = modules.clusters[modules.order[module_idx_to_remove]]

    print("Will remove module index {} (ID: {})".format(
        args.module_remove_idx,
        module_to_remove.id
    ))

    df_entries = []
    for gibbs_idx, alpha, delta, deviation in simulate_random_perturbations(
            study,
            mcmc,
            module_to_remove,
            args.sim_max,
            args.dt
    ):
        df_entries.append({
            'GibbsIdx': gibbs_idx,
            'Alpha': alpha,
            'Delta': delta,
            'Deviation': deviation
        })

    out_path = Path(args.out_path)
    out_path.parent.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(df_entries).to_csv(out_path, index=False, sep='\t')


def simulate_random_perturbations(
        study: md2.Study,
        mcmc: md2.BaseMCMC,
        module: _Cluster,
        sim_max: float,
        dt: float
) -> Iterator[Tuple[int, float, float, float]]:
    alphas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    deltas = [-0.5, -1.0, -1.5, -2.0]
    num_samples = 100

    total_samples = mcmc.n_samples - mcmc.burnin
    stride = total_samples // num_samples

    M = study.matrix(dtype='abs', agg='mean', times='intersection', qpcr_unnormalize=True)
    day21_levels = M[:, 19]

    growths = mcmc.graph[STRNAMES.GROWTH_VALUE].get_trace_from_disk(section="posterior")
    self_interactions = mcmc.graph[STRNAMES.SELF_INTERACTION_VALUE].get_trace_from_disk(section="posterior")
    interactions = mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk(section="posterior")
    interactions[np.isnan(interactions)] = 0
    self_interactions = -np.absolute(self_interactions)
    for i in range(self_interactions.shape[1]):
        interactions[:, i, i] = self_interactions[:, i]

    for gibbs_idx, alpha, delta in itertools.product(
            range(0, total_samples, stride),
            alphas,
            deltas
    ):
        init, r, A = exclude_cluster_from(
            module,
            study.taxa,
            day21_levels,
            growths[gibbs_idx],
            interactions[gibbs_idx]
        )
        x_baseline = run_fwsim_no_pert(
            growth=r, interactions=A,
            initial=init[:, None], sim_max=sim_max, dt=dt, n_days=64
        )

        perts = apply_random_perts(r.shape[0], alpha, delta)
        x_pert = run_fwsim(
            growth=r, interactions=A, pert_strengths=perts,
            pert_start=21, pert_end=34,
            initial=init[:, None],
            sim_max=sim_max, dt=dt, n_days=64
        )

        yield gibbs_idx, alphas, deltas, compute_deviation(x_baseline, x_pert)


def compute_deviation(x1: np.ndarray, x2: np.ndarray, eps: float = 1e5) -> float:
    return float(np.mean(np.log10(x1 + eps) - np.log10(x2 + eps)))


def apply_random_perts(n_taxa: int, fraction: float, strength: float) -> np.ndarray:
    perts = np.zeros(n_taxa, dtype=float)
    n_otus_perturb = max(0, int(n_taxa * fraction))
    selection = np.random.choice(a=n_taxa, size=n_otus_perturb, replace=False)
    perts[selection] = strength
    return perts


def run_fwsim(growth, interactions, pert_strengths, pert_start, pert_end, initial, sim_max, dt, n_days):
    dyn = md2.model.gLVDynamicsSingleClustering(
        growth=growth,
        interactions=interactions,
        perturbations=pert_strengths,
        perturbation_starts=[pert_start],
        perturbation_ends=[pert_end],
        start_day=0,
        sim_max=sim_max
    )

    x = md2.integrate(
        dynamics=dyn,
        initial_conditions=initial,
        dt=dt,
        n_days=n_days,
        subsample=False
    )
    fwsim_values = x['X']
    return fwsim_values


def run_fwsim_no_pert(growth, interactions, initial, sim_max, dt, n_days):
    dyn = md2.model.gLVDynamicsSingleClustering(
        growth=growth,
        interactions=interactions,
        perturbations=None,
        perturbation_starts=[],
        perturbation_ends=[],
        start_day=0,
        sim_max=sim_max
    )

    x = md2.integrate(
        dynamics=dyn,
        initial_conditions=initial,
        dt=dt,
        n_days=n_days,
        subsample=False
    )
    fwsim_values = x['X']
    return fwsim_values


def exclude_cluster_from(
        cluster_to_remove: _Cluster,
        taxa: md2.TaxaSet,
        initial_conditions: np.ndarray,
        growth_rates: np.ndarray,
        interactions: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    oidx_to_remove = set(oidx for oidx in cluster_to_remove.members)
    oidx_to_keep = [
        oidx
        for oidx in range(len(taxa))
        if oidx not in oidx_to_remove
    ]

    return (
        initial_conditions[oidx_to_keep],
        growth_rates[oidx_to_keep],
        interactions[np.ix_(oidx_to_keep, oidx_to_keep)]
    )


if __name__ == "__main__":
    main()
