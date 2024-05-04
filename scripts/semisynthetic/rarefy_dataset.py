from pathlib import Path

import argparse
import numpy as np
import mdsine2 as md2


def rarify_counts(counts: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    if counts.ndim != 1:
        raise ValueError("Only 1-D count vectors are supported.")

    starting_total = np.sum(counts)
    if starting_total < n:
        print("Cannot rarify a counts vector of total = {} into a target total of {}. Applying a ceiling of {}.".format(
            starting_total, n, starting_total
        ))
        n = starting_total

    if starting_total == n:
        return counts
    else:
        nz = counts.nonzero()[0]
        unpacked = np.concatenate([np.repeat(np.array(i, ), counts[i]) for i in nz])
        permuted = rng.permutation(unpacked)[:n]

        result = np.zeros(len(counts), dtype=int)
        for p in permuted:
            result[p] += 1

        return result


def main(
        study_path: Path,
        out_path: Path,
        read_depth: int,
        seed: int
):
    """
    Rarify a dataset pickle file.
    """
    source_study = md2.Study.load(str(study_path))
    new_study = md2.Study(taxa=source_study.taxa)
    rng = np.random.default_rng(seed)

    for source_subj in source_study:
        new_study.add_subject(name=source_subj.name)
        new_subj = new_study[source_subj.name]
        for tidx, t in enumerate(source_subj.times):
            new_subj.reads[t] = rarify_counts(source_subj.reads[t], read_depth, rng)
            new_subj.qpcr[t] = source_subj.qpcr[t]
        new_subj.times = source_subj.times
    new_study.perturbations = source_study.perturbations
    new_study.save(str(out_path))
    print(f"Rarified: {study_path} -> {out_path}")


# ============ CLI interface
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--study', '-s', dest="study_path",
        type=str, required=True,
        help="The path to the original (real) datas Study pickle file."
    )
    parser.add_argument(
        '--out', '-o', dest='out_path',
        type=str, required=True,
        help="The desired full path to which the synthetic study object will be saved."
    )
    parser.add_argument(
        '--read-depth', '-r', dest='read_depth',
        type=int, required=True,
        help="The overall read depth to simulate per timepoint."
    )
    parser.add_argument(
        '--seed', dest='seed',
        type=int, required=True,
        help="The random seed to use for sampling randomness."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        study_path=Path(args.study_path),
        out_path=Path(args.out_path),
        read_depth=args.read_depth,
        seed=args.seed,
    )
