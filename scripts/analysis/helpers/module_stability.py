import argparse
import itertools
from pathlib import Path
from typing import Tuple, Iterator, Optional, Set

import numpy as np
import pandas as pd

import mdsine2 as md2
from mdsine2 import Clustering, generate_cluster_assignments_posthoc
from mdsine2.names import STRNAMES

from tqdm import tqdm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mcmc-path', '-m', type=str, dest='mcmc_path', required=True)
    parser.add_argument('--study', '-s', dest='study', type=str, required=True,
                        help="The path to the relevant Study object containing the input data (subjects, taxa).")
    parser.add_argument('--module-remove-idx', '-i', dest='module_remove_idx', type=int, required=False,
                        help='Specify which module to remove, specified by index (zero-indexed)')
    parser.add_argument('--n_module_replicates', '-n', type=int, required=True,
                        help='Specify the number of replicate modules to use.')

    parser.add_argument('--out-path', '-o', dest='out_path', required=True, type=str)

    parser.add_argument('--seed', required=False, type=int, default=31415)
    parser.add_argument('--num-trials', dest='n_trials', required=False, type=int, default=100)
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

    if module_idx_to_remove is None:
        module_to_remove = None
    else:
        modules: Clustering = load_modal_clustering(mcmc)
        module_to_remove = modules.clusters[modules.order[module_idx_to_remove]]
        print("Will remove module index {} (ID: {})".format(
            args.module_remove_idx,
            module_to_remove.id
        ))

    df_entries = []
    print("Computing sims for module.")
    for gibbs_idx, alpha, delta, trial, deviation in simulate_random_perturbations(
            study,
            mcmc,
            set(oidx for oidx in module_to_remove.members),
            args.sim_max,
            args.dt,
            args.n_trials
    ):
        n_perturbed = int(len(study.taxa) * alpha)
        df_entries.append({
            'ModuleType': 'Fixed',
            'ModuleId': f'F{module_idx_to_remove}',
            'GibbsIdx': gibbs_idx,
            'Alpha': alpha,
            'Delta': delta,
            'Trial': trial,
            'NumPerturbed': n_perturbed,
            'Deviation': deviation
        })

    print("Computing sims for random replicates.")
    for replicate_idx in tqdm(range(args.n_module_replicates)):
        module_random = np.random.choice(
            a=len(study.taxa),
            size=len(module_to_remove.members),
            replace=False
        )
        for gibbs_idx, alpha, delta, trial, deviation in simulate_random_perturbations(
                study,
                mcmc,
                module_random,
                args.sim_max,
                args.dt,
                args.n_trials
        ):
            n_perturbed = int(len(study.taxa) * alpha)
            df_entries.append({
                'ModuleType': 'Random',
                'ModuleId': f'R{replicate_idx}',
                'GibbsIdx': gibbs_idx,
                'Alpha': alpha,
                'Delta': delta,
                'Trial': trial,
                'NumPerturbed': n_perturbed,
                'Deviation': deviation
            })

    out_path = Path(args.out_path)
    out_path.parent.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(df_entries).to_csv(out_path, index=False, sep='\t')


def load_modal_clustering(mcmc) -> Clustering:
    clustering = mcmc.graph[STRNAMES.CLUSTERING_OBJ]
    ret = generate_cluster_assignments_posthoc(clustering, n_clusters='mode', set_as_value=False)
    clustering.from_array(ret)
    return clustering


def simulate_random_perturbations(
        study: md2.Study,
        mcmc: md2.BaseMCMC,
        module: Optional[Set[int]],
        sim_max: float,
        dt: float,
        n_trials: int
) -> Iterator[Tuple[int, float, float, int, float]]:
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

    trials = list(range(1, n_trials + 1, 1))
    for gibbs_idx in tqdm(range(0, total_samples, stride), total=(total_samples // stride)):
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

        for alpha, delta, trial in itertools.product(
                alphas,
                deltas,
                trials
        ):
            perts = apply_random_perts(r.shape[0], alpha, delta)
            x_pert = run_fwsim(
                growth=r, interactions=A, pert_strengths=perts,
                pert_start=21, pert_end=34,
                initial=init[:, None],
                sim_max=sim_max, dt=dt, n_days=64
            )

            yield gibbs_idx, alpha, delta, trial, compute_deviation(x_baseline, x_pert, dt=dt)


def compute_deviation(x1_traj: np.ndarray, x2_traj: np.ndarray, dt: float, eps: float = 1e5) -> float:
    n = int(0.5 / dt)  # number of timepoints to average over.
    x1 = x1_traj[:, -n:].mean(axis=1)
    x2 = x2_traj[:, -n:].mean(axis=1)
    return float(
        np.mean(
            np.abs(
                np.log10(x1 + eps) - np.log10(x2 + eps)
            )
        )
    )


def apply_random_perts(n_taxa: int, fraction: float, strength: float) -> np.ndarray:
    """
    Apply the specified perturbation to a random fraction of taxa.
    :param n_taxa:
    :param fraction:
    :param strength:
    :return:
    """
    perts = np.zeros(n_taxa, dtype=float)
    n_otus_perturb = max(0, int(n_taxa * fraction))
    selection = np.random.choice(a=n_taxa, size=n_otus_perturb, replace=False)
    perts[selection] = strength
    return perts


def run_fwsim(growth, interactions, pert_strengths, pert_start, pert_end, initial, sim_max, dt, n_days):
    dyn = md2.model.gLVDynamicsSingleClustering(
        growth=growth,
        interactions=interactions,
        perturbations=[pert_strengths],
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
        cluster_to_remove: Optional[Set[int]],
        taxa: md2.TaxaSet,
        initial_conditions: np.ndarray,
        growth_rates: np.ndarray,
        interactions: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if cluster_to_remove is None:
        return initial_conditions, growth_rates, interactions

    oidx_to_keep = [
        oidx
        for oidx in range(len(taxa))
        if oidx not in cluster_to_remove
    ]

    return (
        initial_conditions[oidx_to_keep],
        growth_rates[oidx_to_keep],
        interactions[np.ix_(oidx_to_keep, oidx_to_keep)]
    )


if __name__ == "__main__":
    main()
