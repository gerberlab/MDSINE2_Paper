import argparse
import numpy as np
import mdsine2 as md2
from pathlib import Path


def generate_initial_conditions(
        real_data_study: md2.Study,
        n_subjects: int,
        rng: np.random.Generator
) -> np.ndarray:
    """
    Sample initial conditions for N subjects across T taxa.
    :param real_data_study:
    :param n_subjects:
    :param rng:
    :return: A (N x T) array of initial condition values. This will be used downstream to forward-simulate.
    """
    x_all = []
    for subj in real_data_study:
        x = subj.matrix()['abs'][:, 0]  # the first timeslice.
        x_all.append(x)
    x_all = np.stack(x_all, axis=0)  # (n_real_mice x n_taxa)

    # fit lognormal distribution for each taxa.
    x_log = np.log(x_all)
    sigma_log = np.std(x_log, axis=0)
    n_taxa = len(real_data_study.taxa)

    sampled_abunds = x_all[0] * rng.lognormal(mean=np.zeros(n_taxa), sigma=sigma_log, size=(n_subjects, n_taxa))
    sampled_abunds[np.isnan(sampled_abunds)] = 0.0
    return sampled_abunds + 1e5


def main(real_data_study: md2.Study, n_subjects: int, seed: int, out_path: Path):
    rng: np.random.Generator = np.random.default_rng(seed)
    initial_conditions = generate_initial_conditions(real_data_study, n_subjects, rng)
    np.save(out_path, initial_conditions)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--real-data-pkl", dest="real_data_pkl", type=str, required=True)
    parser.add_argument("-n", "--n-subjects", dest="n_subjects", type=int, required=True)
    parser.add_argument("-s", "--seed", dest="seed", type=int, required=True)
    parser.add_argument("-o", "--outpath", dest="outpath", type=str, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        real_data_study=md2.Study.load(args.real_data_pkl),
        n_subjects=args.n_subjects,
        seed=args.seed,
        out_path=Path(args.outpath)
    )