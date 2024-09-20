import argparse
import mdsine2 as md2
from pathlib import Path
from typing import Tuple, List
import numpy as np
import pickle as pkl
from base import GLVParamSet, forward_simulate_subject
import pandas as pd


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


def forward_simulate(
        glv_params: GLVParamSet,
        pert_starts: List[float],
        pert_ends: List[float],
        n_subjects: int,
        n_perts: int,
        rng: np.random.Generator,
        end_t: float, sim_dt: float,
        sim_max: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit a log-normal distribution using real data. Then uses this distribution to sample initial conditions using the "rng" object.
    For each subject (1 through n_subjects), a separate initial condition is sampled iid.
    Then forward-simulation is performed for each timepoint (0, dt, 2*dt, ..., end_t) using the GLV params from the MDSINE2 code.
    """
    # TODO implement this.
    raise NotImplementedError("TODO -- implement this function!")

    # BEGIN starter code below. Use this, or start over somehow.
    n_taxa = len(glv_params.growth)
    trajectories = []
    initial_conditions: np.ndarray = rng.random(size=(n_subjects, n_taxa))  # TODO: need to sample initial conditions using "rng". Maybe need to fit a log-normal distribution.
    for subj in range(n_subjects):
        x, t = forward_simulate_subject(
            glv_params, pert_starts, pert_ends,
            initial_conditions=initial_conditions,
            dt=sim_dt,
            sim_max=sim_max,
            time_points=[1., end_t],  # BE CAREFUL HERE: should simulation start at 1.0 or 0.0?
            time_subset=False  # this has no effect in this scenario. Just keep "false".
        )
    trajectories = np.stack(trajectories, axis=0)
    return trajectories, t
    # END starter code


def truncate_perturbation_list(
        real_data_study: md2.Study,
        glv_params: GLVParamSet,
        n_target_perts: int
) -> Tuple[List[float], List[float], GLVParamSet]:
    """
    Truncate the pert strengths from glv_params and also read off the timepoints for the perturbations.
    :param real_data_study:
    :param glv_params:
    :param n_target_perts:
    :return:
    """
    if len(n_target_perts) < 0:
        raise ValueError("n_target_perts must be a nonnegative number.")
    elif len(n_target_perts) > glv_params.perturbations.shape[0]:  # TODO: Check that perturbation strength matrix has shape (n_perts, n_taxa)
        raise ValueError("n_target_perts cannot exceed the number of existing perturbations from real data ({})".format(
            glv_params.perturbations.shape[0]
        ))

    pert_starts: List[float] = ??  # TODO: create a new pert_starts array according to n_perts
    pert_ends: List[float] = ??  # TODO: create a new pert_ends array according to n_perts
    glv_params: GLVParamSet = ??  # TODO: create a new glv_params object according to n_perts
    return pert_starts, pert_ends, glv_params


def generate_perturbation_table(
        glv_params: GLVParamSet,
        pert_starts: List[float],
        pert_ends: List[float],
        n_subjects: int
) -> pd.DataFrame:
    """
    Create a dataframe in the same format as the real data "perturbations.tsv" file.
    :return:
    """
    raise NotImplementedError("TODO!")

    # Example code
    df_entries = []
    df_entries.append({
        "name": "dummy_pert_name",
        "start": 0.0,
        "end": 0.0,
        "subject": f"dummy_mouse_name_{n_subjects}"
    })
    return pd.DataFrame(df_entries)


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

    pert_starts, pert_ends, glv_params = truncate_perturbation_list(real_data_study, glv_params, n_perts)
    pert_df = generate_perturbation_table(glv_params, pert_starts, pert_ends, n_subjects)
    trajectories, full_dense_timepoints = forward_simulate(real_data_study, glv_params, n_subjects, n_perts, rng, end_t, sim_dt, sim_max)

    assert full_dense_timepoints[-1] == end_t
    outdir.mkdir(exist_ok=True, parents=True)
    np.save(outdir / "trajectories.npy", trajectories)
    np.save(outdir / "timepoints.npy", full_dense_timepoints)
    pert_df.to_csv(outdir / "perturbations.tsv", index=False)


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
        sim_max=args.sim_max
    )
