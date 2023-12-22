from pathlib import Path
from typing import Tuple

import argparse
import numpy as np

import mdsine2 as md2
from mdsine2.names import STRNAMES


# ======= Helper functions
def create_synthetic_dataset(
        source_study: md2.Study,
        fixed_cluster_mcmc: md2.BaseMCMC,
        sim_read_depth: int,
        dirichlet_alpha: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator
) -> md2.Study:
    synthetic_study = md2.Study(taxa=source_study.taxa)
    synthetic_study.perturbations = source_study.perturbations
    read_samples_per_subj, qpcr_samples_per_subj = sample_data_from_posterior_trajectory(
        source_study,
        fixed_cluster_mcmc,
        sim_read_depth,
        dirichlet_alpha,
        qpcr_noise_scale,
        rng=rng
    )

    for source_subj in source_study:
        new_subj = synthetic_study.add_subject(name=source_subj.name)

        # Add timepoints
        new_subj.times = source_subj.times

        # Add reads and biomass
        for tidx, t in enumerate(new_subj.times):
            sampled_reads = read_samples_per_subj[source_subj.name]
            sampled_qpcr = qpcr_samples_per_subj[source_subj.name]
            new_subj.reads[t] = sampled_reads[:, tidx]
            new_subj.qpcr[t] = md2.qPCRdata(cfus=sampled_qpcr[tidx], mass=1., dilution_factor=1.)
    return synthetic_study


def sample_data_from_posterior_trajectory(
        study: md2.Study,
        mcmc: md2.BaseMCMC,
        read_depth: int,
        dirichlet_alpha: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator
) -> Tuple[np.ndarray, np.ndarray]:
    """
    :return: A tuple.
    (1) A dictionary of (Taxa x Timepoint) abundance samples per subject,
    (2) A dictionary of (Timepoint x Replicate) qPCR samples per subject.
    """
    trajectory_set: md2.posterior.TrajectorySet = mcmc.graph[STRNAMES.LATENT_TRAJECTORY]
    reads_per_subject = {}
    qpcr_per_subject = {}
    for subj_idx, subj in enumerate(study):
        traj_posterior_samples = trajectory_set.value[subj_idx].get_trace_from_disk(section='posterior')
        n_samples, n_taxa, n_union_timepoints = traj_posterior_samples.shape

        if n_taxa != len(study.taxa):
            raise ValueError("Matrix's n_taxa dimension ({}) doesn't match taxaset size ({})".format(
                n_taxa, len(study.taxa)
            ))

        if n_union_timepoints != len(subj.times):
            print("Matrix's n_timepoints dimension ({}) doesn't match subject `{}` timepoints ({}). Reverse timepoint lookup is required.".format(
                n_union_timepoints, subj.name, len(subj.times)
            ))

            """
            Dec-22-2023
            
            According to code by @dkaplan, MDSINE2 creates "fake" timepoints via linear interpolation for 
            some of the subjects, so the matrix shape might be bigger than the data's provided timepoint array shape.
            
            Refer to the "essential_timepoints" setting of the "Filtering" class, whose default implementation is to
            set each subject's timepoint collection to the union of all the subjects. 
            """
            _union_indices = {t: i for i, t in enumerate(mcmc.graph.data.times[subj_idx])}  # indices of the LARGE timepoint array (with fake timepoints)
            _t_indices = [_union_indices[t] for t in subj.times]  # The ORIGINAL timepoints from the subject's raw measurements
            traj_posterior_samples = traj_posterior_samples[:, :, _t_indices]  # the submatrix containing only the target measurement timepoints.

        read_counts = np.array([
            sample_reads(read_depth, trajs, dirichlet_alpha, rng=rng)
            for trajs in traj_posterior_samples
        ]) # n_samples x n_taxa x n_timepoints
        reads_per_subject[subj.name] = np.floor(np.mean(read_counts, axis=0))  # reads should be integers.

        qpcrs = np.array([
            sample_qpcr(np.sum(trajs, axis=0), qpcr_noise_scale, n_replicates=3, rng=rng)
            for trajs in traj_posterior_samples
        ])
        qpcr_per_subject[subj.name] = np.mean(qpcrs, axis=0)
    return reads_per_subject, qpcr_per_subject


def sample_reads(
        total_reads: int,
        abundance_trajectories: np.ndarray,
        dirichlet_alpha: float,
        rng: np.random.Generator
) -> np.ndarray:
    """
    Sample reads from the (Taxa, Timepoints)-shaped array of trajectories.
    :param total_reads: The read depth to simulate.
    :param abundance_trajectories: The (Taxa, Timepoints)-shaped array of trajectories.
    :param dirichlet_alpha: The dirichlet "alpha" scaling parameter, defined as \alpha := \sum_i \alpha_i.
    :param rng: The numpy Generator object to use (for reproducibility with seeds)
    :return: Abundances in the shape (Taxa, Timepoints).
    """
    _, n_timepoints = abundance_trajectories.shape
    total_mass = np.sum(abundance_trajectories, axis=0)
    reads_all = []  # A list of (Taxa,)-shaped arrays.
    for tidx in range(n_timepoints):
        rel_abund = abundance_trajectories[:, tidx] / total_mass[tidx]
        alpha = dirichlet_alpha * rel_abund
        reads_all.append(dirichlet_multinomial(alpha, total_reads, rng))
    return np.stack(reads_all, axis=1)


def dirichlet_multinomial(alpha: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    alpha_positive = alpha[alpha > 0]
    x = np.zeros(len(alpha), dtype=int)
    x[alpha > 0] = rng.multinomial(n, rng.dirichlet(alpha_positive))
    return x


def sample_qpcr(
        total_mass: np.ndarray,
        noise_scale: float,
        n_replicates: int,
        rng: np.random.Generator
) -> np.ndarray:
    """
    Sample qPCR measurements for the specified biomass trajectory.
    :param total_mass: A 1-D array specifying biomass for each timepoint.
    :param noise_scale: The log-normal noise scale (std-deviation) parameter.
    :param n_replicates: The number of replicates (e.g. 3 for triplicates) of qpcr measurements.
    :param rng: The numpy Generator object to use (for reproducibility with seeds)
    """
    n_times, = total_mass.shape
    return np.exp(
        np.expand_dims(np.log(total_mass), axis=1)  # shape (T, 1)
        +
        noise_scale * rng.normal(loc=0.0, scale=1.0, size=(n_times, n_replicates))
    )


# ============ CLI interface
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--study', '-s', dest="study_path",
        type=str, required=True,
        help="The path to the original (real) datas Study pickle file."
    )
    parser.add_argument(
        '--fixed-module-pkl', '-fm', dest='fixed_module_pkl_path',
        type=str, required=True,
        help="The path to the pickled MCMC run of the fixed-module inference."
    )
    parser.add_argument(
        '--read-depth', '-r', dest='read_depth',
        type=int, required=True,
        help="The overall read depth to simulate per timepoint."
    )
    parser.add_argument(
        '--alpha', '-a', dest='dirichlet_alpha',
        type=float, required=True,
        help="The alpha parameter of the dirichlet distribution (the overall scale: sum_i alpha_i)"
    )
    parser.add_argument(
        '--qpcr-noise-scale', '-q', dest='qpcr_noise_scale',
        type=float, required=True,
        help="The qPCR noise scale (geometric noise stdev scaling)"
    )
    parser.add_argument(
        '--seed', dest='seed',
        type=int, required=True,
        help="The random seed to use for sampling randomness."
    )
    return parser.parse_args()


def main(
        study_path: Path,
        fixed_module_pkl_path: Path,
        read_depth: int,
        dirichlet_alpha: int,
        qpcr_noise_scale: float,
        seed: int
):
    """
    Read the fixed-module inference, use the stored "filtered-state" (latent traj X) samples.

    :param fixed_module_pkl_path:
    :param read_depth:
    :return:
    """
    rng = np.random.default_rng(seed)
    create_synthetic_dataset(
        source_study=md2.Study.load(str(study_path)),
        fixed_cluster_mcmc=md2.BaseMCMC.load(str(fixed_module_pkl_path)),
        sim_read_depth=read_depth,
        dirichlet_alpha=dirichlet_alpha,
        qpcr_noise_scale=qpcr_noise_scale,
        rng=rng
    )


if __name__ == "__main__":
    args = parse_args()
    main(
        study_path=Path(args.study_path),
        fixed_module_pkl_path=Path(args.fixed_module_pkl_path),
        read_depth=args.read_depth,
        dirichlet_alpha=args.dirichlet_alpha,
        qpcr_noise_scale=args.qpcr_noise_scale,
        seed=args.seed
    )
