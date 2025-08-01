"""
generates the data in a format that is compatible to work with the original cLV code
"""
import numpy as np
import pickle as pkl
import argparse
from pathlib import Path

import mdsine2 as md2


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--mdsine2_study", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    sample_id_to_subject_id = {}
    subject_id_time = {}
    subject_id_u = {}

    args = parse_arguments()
    study = md2.Study.load(args.mdsine2_study)
    taxa = study.taxa

    subject_counts = {
        subj.name: subj.matrix()["abs"] for subj in study
    }

    Y = []
    U = []
    T = []
    zero_counts = 0
    total_counts = 0

    for subj in study:
        counts = subj.matrix()["abs"].T
        times = subj.times
        perturb_ids = np.zeros(len(times), dtype=int)
        if study.perturbations is not None:
            for p_idx, pert in enumerate(study.perturbations):
                start = pert.starts[subj.name]
                end = pert.ends[subj.name]
                indices, = np.where((times >= start) & (times <= end))
                perturb_ids[indices] = p_idx + 1

        assert len(times) == counts.shape[0]
        assert counts.shape[1] == len(taxa)

        zero_counts += np.sum(counts == 0)
        total_counts += counts.size
        Y.append(counts)
        U.append(perturb_ids.reshape((-1, 1)))
        T.append(times)

    out_dir = Path(args.output_dir)
    with open(out_dir / "Y.pkl", "wb") as f:
        pkl.dump(Y, f)
    with open(out_dir / "U.pkl", "wb") as f:
        pkl.dump(U, f)
    with open(out_dir / "T.pkl", "wb") as f:
        pkl.dump(T, f)
