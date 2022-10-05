import argparse
import itertools
from pathlib import Path
from typing import Tuple, Iterator, Optional, Set, List

import numpy as np
import scipy.stats
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

import mdsine2 as md2
from tqdm import tqdm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs-dir', '-i', type=str, dest='inputs_dir', required=True)
    parser.add_argument('--study', '-s', dest='study', type=str, required=True,
                        help="The path to the relevant Study object containing the input data (subjects, taxa).")
    parser.add_argument('--module-remove-idx', '-m', dest='module_remove_idx', type=int, required=True,
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
    inputs_dir = Path(args.inputs_dir)
    module_idx_to_remove = args.module_remove_idx

    modules: List[List[int]] = load_modal_clustering(inputs_dir)
    module_to_remove = modules[module_idx_to_remove]
    print("Will remove module index {} (Size {})".format(
        args.module_remove_idx,
        len(module_to_remove)
    ))

    df_entries = []
    print("Computing sims for module.")
    for alpha, delta, trial, gibbs_idx, deviation in simulate_random_perturbations(
            study,
            inputs_dir,
            set(oidx for oidx in module_to_remove),
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
            size=len(module_to_remove),
            replace=False
        )
        for alpha, delta, trial, gibbs_idx, deviation in simulate_random_perturbations(
                study,
                inputs_dir,
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


def load_modal_clustering(inputs_dir: Path) -> List[List[int]]:
    agglomeration_file = inputs_dir / "agglomeration.npy"
    if agglomeration_file.exists():
        agglom = np.load(str(agglomeration_file))
    else:
        A = 1 - np.load(str(inputs_dir / "coclusters.npy"))
        n = scipy.stats.mode(np.load(str(inputs_dir / "n_clusters.npy")))[0][0]

        linkage = 'average'
        c = AgglomerativeClustering(
            n_clusters=n,
            affinity='precomputed',
            linkage=linkage
        )

        agglom = c.fit_predict(A)
        np.save(str(agglomeration_file), agglom)

    clusters = []
    for cidx in range(np.max(agglom) + 1):  # Make sure to do the (+1) to count the last module.
        cluster = list(np.where(agglom == cidx)[0])
        clusters.append(cluster)
    return clusters


def load_parameters(inputs_dir: Path):
    growths = np.load(str(inputs_dir / "growth.npy"))
    interactions = np.load(str(inputs_dir / "interactions.npy"))
    return growths, interactions


def simulate_random_perturbations(
        study: md2.Study,
        inputs_dir: Path,
        module: Optional[Set[int]],
        sim_max: float,
        dt: float,
        n_trials: int
) -> Iterator[Tuple[float, float, int, int, float]]:
    alphas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    deltas = [-0.5, -1.0, -1.5, -2.0]

    growths, interactions = load_parameters(inputs_dir)
    total_samples = growths.shape[0]

    M = study.matrix(dtype='abs', agg='mean', times='intersection', qpcr_unnormalize=True)
    day21_levels = M[:, 19]

    trials = list(range(1, n_trials + 1, 1))
    for trial in trials:
        gibbs_idx = np.random.choice(total_samples)
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

        for alpha, delta in itertools.product(
                alphas,
                deltas
        ):
            perts = apply_random_perts(r.shape[0], alpha, delta)
            x_pert = run_fwsim(
                growth=r, interactions=A, pert_strengths=perts,
                pert_start=21, pert_end=34,
                initial=init[:, None],
                sim_max=sim_max, dt=dt, n_days=64
            )

            yield alpha, delta, trial, gibbs_idx, compute_deviation(x_baseline, x_pert, dt=dt)


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
