import argparse
from pathlib import Path


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
        fwsim_dir=Path(args.fwsim_dir)
    )