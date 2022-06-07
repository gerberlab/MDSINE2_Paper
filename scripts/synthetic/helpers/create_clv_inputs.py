"""
generates the data in a format that is compatible to work with the original cLV code
"""
import numpy as np
import pickle as pkl
import argparse
from pathlib import Path

import mdsine2 as md2
from scipy.interpolate import UnivariateSpline


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--mdsine2_study", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    return parser.parse_args()


def compute_breakpoints(effects):
    """Break the spline when external perturbations occur."""
    breakpoints = []
    for u in effects:
        in_perturb = False
        v = []
        for i,ut in enumerate(u):
            if i == 0 or i == u.shape[0]-1:
                v.append(1)
                continue
            if np.any(ut) > 0 and not in_perturb:
                v.append(1)
                in_perturb = True
            elif np.any(ut) > 0 and in_perturb:
                v.append(0)
            elif np.all(ut) == 0 and in_perturb:
                i = 1 if v[i-1] == 0 else 0
                v.append(i)
                in_perturb = False
            else:
                v.append(0)
        v = np.array(v)
        breakpoints.append(np.nonzero(v)[0])
    return breakpoints


def denoise(counts, t_pts, effects=None):
    """Takes a sequence of counts at t_pts, and returns denoised estimates
    of latent trajectories."""
    ntaxa = counts[0].shape[1]
    denoised_traj = []

    if effects is not None:
        breakpoints = compute_breakpoints(effects)
        for c,t,b in zip(counts,t_pts,breakpoints):
            denoised = np.zeros(c.shape)
            mass = c.sum(axis=1,keepdims=True)
            p = c / c.sum(axis=1,keepdims=True)
            p[p==0] = 1e-5
            p /= p.sum(axis=1,keepdims=True)
            c = (mass.T*p.T).T
            for i in range(ntaxa):
                for j in range(1,b.size):
                    start = b[j-1]
                    end = b[j]+1
                    k = 5 if end - start <= 3 else 5
                    f = UnivariateSpline(t[start:end],c[start:end,i],k=k)
                    denoised[start:end,i] = f(t[start:end])
            denoised[0] = c[0]
            denoised = np.clip(denoised, np.min(denoised[denoised > 0]), np.inf)
            denoised_traj.append(denoised)
    else:
        for c, t in zip(counts,t_pts):
            denoised = np.zeros(c.shape)
            k = 3 if t.shape[0] <= 5 else 5
            for i in range(ntaxa):
                f = UnivariateSpline(t,c[:,i],k=k)
                denoised[:,i] = f(t)
            denoised = np.clip(denoised, np.min(denoised[denoised > 0]), np.inf)
            denoised /= denoised.sum(axis=1,keepdims=True)
            denoised_traj.append(denoised)

    return denoised_traj


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
        perturb_ids = np.zeros(times.shape)

        assert len(times) == counts.shape[0]
        assert counts.shape[1] == len(taxa)

        zero_counts += np.sum(counts == 0)
        total_counts += counts.size
        Y.append(counts)
        U.append(perturb_ids)
        T.append(times)

    Y_denoised = denoise(Y, T, effects=U)
    out_dir = Path(args.output_dir)
    with open(out_dir / "Y.pkl", "wb") as f:
        pkl.dump(Y, f)
    with open(out_dir / "Y_denoised.pkl", "wb") as f:
        pkl.dump(Y_denoised, f)
    with open(out_dir / "U.pkl", "wb") as f:
        pkl.dump(U, f)
    with open(out_dir / "T.pkl", "wb") as f:
        pkl.dump(T, f)
