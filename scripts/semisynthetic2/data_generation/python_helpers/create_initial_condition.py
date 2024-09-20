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

    """ TODO: create initial conditions in a smarter way. """
    # Dummy code: copy from mouse 1, repeat N times. Add 1e4 to pad the zeroes.
    n_taxa = len(real_data_study.taxa)
    target_real_mouse = real_data_study['2']
    return np.stack([
        target_real_mouse.matrix()['abs'][:, 0] + 1e4
        for _ in range(n_subjects)
    ], axis=0)


def main(real_data_study: md2.Study, n_subjects: int, seed: int, out_path: Path):
    rng: np.random.Generator = np.random.default_rng(seed)
    initial_conditions = generate_initial_conditions(real_data_study, n_subjects, rng)
    np.save(out_path, initial_conditions)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--real-data-pkl", dest="real_data_pkl", type=int, required=True)
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