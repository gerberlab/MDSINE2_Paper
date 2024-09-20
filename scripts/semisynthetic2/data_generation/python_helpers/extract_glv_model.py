from typing import *
import numpy as np
import pickle as pkl
import argparse
from pathlib import Path

import mdsine2 as md2
from scripts.semisynthetic2.data_generation.python_helpers.base import GLVParamSet, forward_simulate


def extract_glv_model(
        study: md2.Study,
        growths: np.ndarray,
        interactions: np.ndarray,
        perturbations: List[np.ndarray],
        coclusterings: np.ndarray,
        dt: float,
        sim_max: float,
) -> Tuple[GLVParamSet, np.ndarray, Dict[str, np.ndarray]]:
    """Pick a best forward simulation based on error evaluation metric."""

    subsample_every = 1000
    if subsample_every > 1:
        print("Subsampling every {} MCMC samples. (overall, there are {} samples)".format(subsample_every, interactions.shape[0]))
    params = [
        GLVParamSet(
            growth=growths[i],
            interactions=interactions[i],
            perturbations=[p[i] for p in perturbations]
        )
        for i in range(0, interactions.shape[0], subsample_every)
    ]
    coclusterings = coclusterings[::subsample_every]

    # Extract the best parameter set using forward simulations.
    param_fwsim_errors = []
    param_num_surviving_modules = []
    # for param_set, coclust_mat in tqdm(zip(params, coclusterings), total=len(params)):
    iter_idx = 0
    for param_set, coclust_mat in zip(params, coclusterings):
        err, surviving_taxa, surviving_modules, total_modules = evaluate_parameter_fwsim(param_set, study, dt, sim_max, coclust_mat)
        # param_fwsim_errors.append(err)
        if len(surviving_taxa) < len(study.taxa):
            param_fwsim_errors.append(np.inf)
        else:
            param_fwsim_errors.append(err)
        param_num_surviving_modules.append(len(surviving_modules))
        print("[iter={}] Modules survived: {} / {}".format(iter_idx, len(surviving_modules), total_modules))
        iter_idx += subsample_every

    best_idx = np.argmin(param_fwsim_errors)
    if param_fwsim_errors[best_idx] == np.inf:
        raise ValueError("Couldn't find best sim with all surviving taxa")
    print("Choice: index {} (total {} entries)".format(
        best_idx,
        len(params)
    ))
    return params[best_idx], coclusterings[best_idx], forward_simulate(params[best_idx], study, dt, sim_max)


# =========== helpers
def extract_clustering_from_matrix(_mat) -> np.ndarray:
    n_items = _mat.shape[0]
    items_left = set(range(n_items))

    clustering_assignments = np.zeros(n_items, dtype=int)

    c_idx = -1
    while len(items_left) > 0:
        c_idx += 1  # new cluster
        x = items_left.pop()  # member
        clustering_assignments[x] = c_idx

        # get all guys in the same cluster
        to_remove = set()
        for y in items_left:
            if _mat[x,y]:
                clustering_assignments[y] = c_idx
                to_remove.add(y)
        items_left = items_left.difference(to_remove)
    return clustering_assignments


def evaluate_parameter_fwsim(
        params: GLVParamSet,
        study: md2.Study,
        dt: float,
        sim_max: float,
        coclust_matrix: np.ndarray
) -> Tuple[float, Set[int], Set[int], int]:
    """
    Evaluate the fwsim error for each subject (And sum them).
    """
    fwsims = forward_simulate(params, study, dt, sim_max)
    errors = []
    eps = 1e-5
    for subj_name, fwsim_trajs in fwsims.items():
        subj = study[subj_name]
        measurements = subj.matrix()['abs']

        # RMS-log10 error, excluding the points where data says zero.
        # fwsim_trajs: [taxa x timepoints] array
        # measurements: [taxa x timepoints] array, with zeroes in some of the entries

        mask = measurements > 0
        diff_logs = np.log10(fwsim_trajs + eps) - np.log10(measurements + eps)
        subj_err = np.sqrt(np.mean(np.square(diff_logs[mask])))
        errors.append(subj_err)

    # What taxa survived? (Only keep taxa which surpasses 1e5 abundance for some timepoint in some synthetic mouse)
    taxa_max_abundance = np.stack([
        fwsim_trajs.max(axis=-1)  # length=n_taxa, after maximizing across timepoints
        for subj_name, fwsim_trajs in fwsims.items()
    ], axis=0).max(axis=0)  # length=n_taxa, max across subjects

    clust_assignments = extract_clustering_from_matrix(coclust_matrix)
    leftover_module_set = set()
    leftover_taxa_set = set()
    for taxa_idx, taxa in enumerate(study.taxa):
        if taxa_max_abundance[taxa_idx] > 1e5:
            leftover_module_set.add(clust_assignments[taxa_idx])
            leftover_taxa_set.add(taxa_idx)
    total_modules = np.max(clust_assignments) + 1

    print("leftover taxa: {}".format(len(leftover_taxa_set)))
    return np.sum(errors), leftover_taxa_set, leftover_module_set, total_modules


# ============ CLI interface
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--study', '-s', dest="study_path",
        type=str, required=True,
        help="The path to the original (real) datas Study pickle file."
    )
    parser.add_argument(
        '--growths', dest='growth_path',
        type=str, required=True, help="The path to the .npy growth rate MCMC samples."
    )
    parser.add_argument(
        '--interactions', dest='interactions_path',
        type=str, required=True, help="The path to the .npy interactions MCMC samples."
    )
    parser.add_argument(
        '--perts', dest='perts_path',
        type=str, required=True, help="The path to the .npz perturbation strengths MCMC samples."
    )
    parser.add_argument(
        '--coclust', dest='coclust_path',
        type=str, required=True, help="The path to the .npz coclustering MCMC samples."
    )
    parser.add_argument(
        '--out-dir', '-o', dest='ground_truth_dir',
        type=str, required=True,
        help="A directory to save the chosen desired ground truth values."
    )
    parser.add_argument('--sim-dt', type=float, required=False, default=0.01)
    parser.add_argument('--sim-max', type=float, required=False, default=1e20)
    return parser.parse_args()


def main(
        study_path: Path,
        growth_path: Path,
        interactions_path: Path,
        perturbations_path: Path,
        coclusterings_path: Path,
        ground_truth_dir: Path,
        sim_max: float,
        sim_dt: float
):
    """
    Evaluate the posterior parameters, and pick the best one.
    :param study_path:
    :param growth_path:
    :param interactions_path:
    :param perturbations_path:
    :param coclusterings_path:
    :param ground_truth_dir:
    :param sim_max:
    :param sim_dt:
    :return:
    """
    source_study = md2.Study.load(str(study_path))
    ground_truth_dir.mkdir(parents=True, exist_ok=True)
    print(f"Using ground truth dir {ground_truth_dir}")

    glv_sim_path = ground_truth_dir / 'glv_best_sim.pkl'

    growths = np.load(growth_path)
    interactions = np.load(interactions_path)
    perturbations_map = np.load(perturbations_path)
    perturbations = [
        perturbations_map[pert.name]
        for pert in source_study.perturbations
    ]
    coclusterings = np.load(coclusterings_path)
    assert coclusterings.shape[0] == interactions.shape[0]  # ensure that the number of samples match

    # Pick the best MCMC sample.
    glv_params, instance_coclusters, fwsims = extract_glv_model(source_study, growths, interactions, perturbations, coclusterings, sim_dt, sim_max)
    with open(glv_sim_path, 'wb') as f:
        pkl.dump(glv_params, f)
    np.save(ground_truth_dir / "coclusters.npy", instance_coclusters)


if __name__ == "__main__":
    args = parse_args()
    main(
        study_path=Path(args.study_path),
        growth_path=Path(args.growth_path),
        interactions_path=Path(args.interactions_path),
        perturbations_path=Path(args.perts_path),
        coclusterings_path=Path(args.coclust_path),
        ground_truth_dir=Path(args.ground_truth_dir),
        sim_max=args.sim_max,
        sim_dt=args.sim_dt
    )

