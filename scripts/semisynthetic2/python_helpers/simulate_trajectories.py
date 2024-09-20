import argparse
import mdsine2 as md2
from pathlib import Path
from typing import Tuple
import numpy as np
import pickle as pkl
from base import GLVParamSet, forward_simulate_subject


def subj_with_most_timepoints(study: md2.Study) -> md2.Subject:
    best_subj = None
    for subj in study:
        if best_subj is None:
            best_subj = subj
        elif len(best_subj.times) < len(subj.times):
            best_subj = subj
    if best_subj is None:
        raise ValueError("No subjects found in study.")
    return best_subj


def forward_simulate(real_data_study: md2.Study, glv_params: GLVParamSet, n_subjects: int, n_perts: int, rng: np.random.Generator, end_t: float, sim_dt: float, sim_max: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit a log-normal distribution using real data. Then uses this distribution to sample initial conditions using the "rng" object.
    For each subject (1 through n_subjects), a separate initial condition is sampled iid.
    Then forward-simulation is performed for each timepoint (0, dt, 2*dt, ..., end_t) using the GLV params from the MDSINE2 code.
    """
    # TODO implement this.
    raise NotImplementedError("TODO -- implement this function!")

    # BEGIN starter code below. Use this, or start over somehow.
    n_taxa = len(real_data_study.taxa)
    initial_conditions = rng.random(size=(n_subjects, n_taxa))
    trajectories = []
    for subj in range(n_subjects):
        x, t = forward_simulate_subject(
            glv_params, pert_starts, pert_ends,  # TODO: need to set the perturbation times correctly.
            initial_conditions,  # TODO: need to sample initial conditions using "rng". Maybe need to fit a log-normal distribution.
            sim_dt,
            sim_max,
            time_points=[1., end_t],  # BE CAREFUL HERE: should simulation start at 1.0 or 0.0?
            time_subset=False  # this has no effect in this scenario. Just keep "false".
        )
    return trajectories, t
    # END starter code


def main(real_data_study: md2.Study, n_subjects: int, n_perts: int, seed: int, outdir: Path, ground_truth_dir: Path, sim_dt: float, sim_max: float):
    """
    Invoke forward_simulate helper function. Save the trajectories/timepoints array into .npy array files on disk.
    """
    target_mouse = subj_with_most_timepoints(real_data_study)
    end_t = target_mouse.times[-1]
    rng = np.random.default_rng(seed=seed)

    # Load ground truth parameters.
    glv_sim_path = ground_truth_dir / 'glv_best_sim.pkl'
    with open(glv_sim_path, 'rb') as f:
        glv_params: GLVParamSet = pkl.load(f)

    trajectories, full_dense_timepoints = forward_simulate(real_data_study, glv_params, n_subjects, n_perts, rng, end_t, sim_dt, sim_max)

    assert full_dense_timepoints[-1] == end_t
    outdir.mkdir(exist_ok=True, parents=True)
    np.save(outdir / "trajectories.npy", trajectories)
    np.save(outdir / "timepoints.npy", full_dense_timepoints)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--real-data-pkl", dest="real_data_pkl", type=int, required=True)
    parser.add_argument("-n", "--n-subjects", dest="n_subjects", type=int, required=True)
    parser.add_argument("-p", "--n-perts", dest="n_perts", type=int, required=True)
    parser.add_argument("-s", "--seed", type=int, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-g", "--ground-truth-dir", dest="ground_truth_dir", type=str, required=True)
    parser.add_argument("-dt", "--sim-dt", dest="sim_dt", type=float, required=False, default=0.01)
    parser.add_argument("-sm", "--sim-max", dest="sim_max", type=float, required=False, default=1e20)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        real_data_study=md2.Study.load(args.real_data_pkl),
        n_subjects=args.n_subjects,
        n_perts=args.n_perts,
        seed=args.seed,
        outdir=Path(args.outdir),
        ground_truth_dir=Path(args.ground_truth_dir),
        sim_dt=args.sim_dt,
    )
