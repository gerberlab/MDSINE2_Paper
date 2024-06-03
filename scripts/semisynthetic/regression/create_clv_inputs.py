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

    maximum_post_rescale_target = 100.0
    n_perts = len(study.perturbations)
    for subj in study:
        abund = subj.matrix()["abs"].T
        abund_rescale_factor = maximum_post_rescale_target / np.max(abund)

        times = subj.times
        perturb_ids = np.zeros((len(times), n_perts), dtype=int)
        for p_idx, pert in enumerate(study.perturbations):
            start = pert.starts[subj.name]
            end = pert.ends[subj.name]
            timepoint_indices, = np.where((times >= start) & (times <= end))
            perturb_ids[timepoint_indices, p_idx] = 1

        assert len(times) == abund.shape[0]
        assert abund.shape[1] == len(taxa)

        zero_counts += np.sum(abund == 0)
        total_counts += abund.size
        Y.append(abund * abund_rescale_factor)
        U.append(perturb_ids)
        T.append(times)

    print("Checker")
    print("len(Y):{}, len(U):{}, len(T):{}".format(len(Y), len(U), len(T)))
    print("Shape Y[0]:{}, Shape U[0]:{}, Shape T[0]:{}".format(Y[0].shape, U[0].shape, T[0].shape))
    print("Rescale factor: {}".format(abund_rescale_factor))

    out_dir = Path(args.output_dir)
    with open(out_dir / "Y.pkl", "wb") as f:
        pkl.dump(Y, f)
    with open(out_dir / "U.pkl", "wb") as f:
        pkl.dump(U, f)
    with open(out_dir / "T.pkl", "wb") as f:
        pkl.dump(T, f)
    with open(out_dir / "rescale.np", "wb") as f:
        np.save(f, np.array(abund_rescale_factor))
