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
        initial_conditions: np.ndarray,
        pert_starts: List[float],
        pert_ends: List[float],
        n_subjects: int,
        end_t: float, sim_dt: float,
        sim_max: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit a log-normal distribution using real data. Then uses this distribution to sample initial conditions using the "rng" object.
    For each subject (1 through n_subjects), a separate initial condition is sampled iid.
    Then forward-simulation is performed for each timepoint (0, dt, 2*dt, ..., end_t) using the GLV params from the MDSINE2 code.
    """
    trajectories = []
    assert n_subjects > 0
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


def truncate_perturbation_list(
        real_mouse_pert_names: List[str],
        real_mouse_timepoints: np.ndarray,
        real_mouse_pert_starts: List[float],
        real_mouse_pert_ends: List[float],
        pert_strengths: np.ndarray,
        n_target_perts: int
) -> Tuple[List[str], List[float], List[float], np.ndarray]:
    """
    Truncate the pert strengths from glv_params and also read off the timepoints for the perturbations.
    :param real_mouse_pert_names: The list of perturbation names. (Call this length "p")
    :param real_mouse_timepoints: The timepoints (an array of timepoint values) for one real mouse of our choice.
    :param real_mouse_pert_starts: The perturbation start times (a list of length p) for one real mouse of our choice.
    :param real_mouse_pert_ends: The perturbation end times (a list of length p) for one real mouse of our choice.
    :param pert_strengths: A (p x n_taxa) array of pert strengths.
    :param n_target_perts: The number of perts to truncate down to.
    :return:
    """
    if n_target_perts < 0:
        raise ValueError("n_target_perts must be a nonnegative number.")
    elif n_target_perts > pert_strengths.shape[0]:  # TODO: Check that perturbation strength matrix has shape (n_perts, n_taxa)
        raise ValueError("n_target_perts cannot exceed the number of existing perturbations from real data ({})".format(
            pert_strengths.shape[0]
        ))

    """ =========================================== TODO -- Task #1 ============================================ """
    pert_names: List[str] = real_mouse_pert_names
    pert_starts: List[float] = real_mouse_pert_starts  # TODO: create a new pert_starts array according to n_perts
    pert_ends: List[float] = real_mouse_pert_ends  # TODO: create a new pert_ends array according to n_perts
    pert_strengths: np.ndarray = pert_strengths  # TODO: create a new glv_params object according to n_perts
    return pert_names, pert_starts, pert_ends, pert_strengths


def generate_perturbation_table(
        pert_names: List[str],
        pert_starts: List[float],
        pert_ends: List[float],
        n_subjects: int
) -> pd.DataFrame:
    """
    Create a dataframe in the same format as the real data "perturbations.tsv" file.
    :param pert_names: A (perts)-length array listing out the names of each perturbation.
    :param perturbation_strengths: A (perts x taxa) shaped array.
    :param pert_starts: A (perts)-length array listing out the start time of each perturbation.
    :param pert_ends: A (perts)-length array listing out the end time of each perturbation.
    :param n_subjects: The number of target subjects. In this simulation, we will have the same perturbation window for each subject.
    :return:
    """

    """ ============================================ TODO -- Task #2 ============================================ """
    # Example code
    df_entries = []
    df_entries.append({
        "name": "dummy_pert_name",
        "start": 0.0,
        "end": 0.0,
        "subject": f"dummy_mouse_name_{n_subjects}"
    })
    return pd.DataFrame(df_entries)


def main(
        real_data_study: md2.Study, n_subjects: int, n_perts: int, initial_conditions: np.ndarray,
        outdir: Path, ground_truth_dir: Path, sim_dt: float, sim_max: float
):
    """
    Invoke forward_simulate helper function. Save the trajectories/timepoints array into .npy array files on disk.
    """
    # Pick our favorite mouse: one with the most timepoints.
    target_mouse = subj_with_most_timepoints(real_data_study)

    # Load ground truth parameters.
    glv_sim_path = ground_truth_dir / 'glv_best_sim.pkl'
    with open(glv_sim_path, 'rb') as f:
        glv_params: GLVParamSet = pkl.load(f)

    # The metadata of the real mouse: perturbation names, pert start times, pert end times.
    real_mouse_end_t = target_mouse.times[-1]
    real_mouse_pert_names = [p.name for p in real_data_study.perturbations]
    real_mouse_pert_starts = [
        pert.starts[target_mouse.name]
        for pert in real_data_study.perturbations
    ]
    real_mouse_pert_ends = [
        pert.ends[target_mouse.name]
        for pert in real_data_study.perturbations
    ]

    # New synthetic information.
    pert_names_for_sim, pert_starts_for_sim, pert_ends_for_sim, pert_strengths_for_sim = truncate_perturbation_list(
        real_mouse_pert_names,
        target_mouse.times,
        real_mouse_pert_starts,
        real_mouse_pert_ends,
        glv_params.perturbations,
        n_perts)
    pert_df = generate_perturbation_table(
        pert_names_for_sim,
        pert_starts_for_sim,
        pert_ends_for_sim,
        n_subjects
    )

    sim_glv_params = GLVParamSet(glv_params.growth, glv_params.interactions, pert_strengths_for_sim)
    trajectories, full_dense_timepoints = forward_simulate(
        sim_glv_params,
        initial_conditions,
        pert_starts_for_sim,
        pert_ends_for_sim,
        n_subjects,
        real_mouse_end_t,
        sim_dt,
        sim_max
    )

    assert full_dense_timepoints[-1] == real_mouse_end_t
    outdir.mkdir(exist_ok=True, parents=True)
    np.save(outdir / "trajectories.npy", trajectories)
    np.save(outdir / "timepoints.npy", full_dense_timepoints)
    pert_df.to_csv(outdir / "perturbations.tsv", index=False, sep='\t')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--real-data-pkl", dest="real_data_pkl", type=int, required=True)
    parser.add_argument("-n", "--n-subjects", dest="n_subjects", type=int, required=True)
    parser.add_argument("-p", "--n-perts", dest="n_perts", type=int, required=True)
    parser.add_argument("-i", "--initial-conditions", dest="init_cond_path", type=str, required=True)
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
        initial_conditions=np.load(args.init_cond_path),
        outdir=Path(args.outdir),
        ground_truth_dir=Path(args.ground_truth_dir),
        sim_dt=args.sim_dt,
        sim_max=args.sim_max
    )
