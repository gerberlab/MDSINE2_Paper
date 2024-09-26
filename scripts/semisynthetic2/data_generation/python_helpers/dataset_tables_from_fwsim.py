from typing import *
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from base import generate_mouse_name

import mdsine2 as md2
from mdsine2 import TaxaSet


# ======= Helper functions
def sample_data_from_fwsim(
        forward_sims: np.ndarray,
        subj_names: List[str],
        sim_timepoints: np.ndarray,
        read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], List[float]]:
    """
    :return: A tuple.
    (1) A dictionary of (Taxa x Timepoint) abundance samples per subject,
    (2) A dictionary of (Timepoint x Replicate) qPCR samples per subject.
    """
    assert len(forward_sims.shape) == 3
    reads_per_subject = {}
    qpcr_per_subject = {}
    n_subjs = forward_sims.shape[0]

    # only take 0.5-discretized timepoints.
    start_t = sim_timepoints[0]
    end_t = sim_timepoints[-1]
    timepoint_subset = set(np.arange(start_t, end_t+0.5, step=0.5))
    timepoint_subset_indices = [i for i, t in enumerate(sim_timepoints) if t in timepoint_subset]

    for subj_idx, subj_name in zip(range(n_subjs), subj_names):
        trajs = forward_sims[subj_idx]
        trajs = trajs[:, timepoint_subset_indices]

        if trajs.shape[1] != len(timepoint_subset):
            raise ValueError("Matrix's n_times dimension ({}) doesn't match desired num of timepoints ({})".format(
                trajs.shape[1], len(timepoint_subset)
            ))

        read_counts = sample_reads(read_depth, trajs, negbin_a0, negbin_a1, rng=rng)  # n_taxa x n_timepoints
        reads_per_subject[subj_name] = read_counts

        qpcrs = sample_qpcr(np.sum(trajs, axis=0), qpcr_noise_scale, n_replicates=3, rng=rng)
        qpcr_per_subject[subj_name] = qpcrs
    return reads_per_subject, qpcr_per_subject, sorted(timepoint_subset)


def sample_replicate_data_from_fwsim(
        forward_sim: np.ndarray,
        sim_timepoints: np.ndarray,
        num_physical_replicates: int,
        read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator,
) -> Tuple[Dict[int, np.ndarray], Dict[int, np.ndarray], List[int]]:
    """
    :return: A tuple.
    (1) A dictionary of (Taxa x Timepoint) abundance samples per REPLICATE (key is replicate index),
    (2) A dictionary of (Timepoint x Replicate) qPCR samples per REPLICATE (key is replicate index).
    (3) A pre-set list of the subset timepoint values. (this function only produces replicate data for integer-valued timepoints.)
    """
    if forward_sim.shape[1] != len(sim_timepoints):
        raise ValueError("Matrix's n_times dimension ({}) doesn't match subj times size ({})".format(
            forward_sim.shape[1], len(sim_timepoints)
        ))

    # Only take integer-valued timepoints.
    start_t = sim_timepoints[0]
    end_t = sim_timepoints[-1]
    timepoint_subset = set(np.arange(start_t, end_t, step=1.))
    timepoint_subset_indices = [i for i, t in enumerate(sim_timepoints) if t in timepoint_subset]

    forward_sim_subset = forward_sim[:, timepoint_subset_indices]
    read_counts_all = {}
    qpcrs_all = {}
    for replicate_idx in range(num_physical_replicates):
        read_counts = sample_reads(read_depth, forward_sim_subset, negbin_a0, negbin_a1, rng=rng)  # n_taxa x n_timepoints
        qpcrs = sample_qpcr(np.sum(forward_sim_subset, axis=0), qpcr_noise_scale, n_replicates=3, rng=rng)
        read_counts_all[replicate_idx] = read_counts
        qpcrs_all[replicate_idx] = qpcrs
    return read_counts_all, qpcrs_all, [int(t) for t in sorted(timepoint_subset)]


def sample_reads(
        total_reads: int,
        abundance_trajectories: np.ndarray,
        negbin_a0: float,
        negbin_a1: float,
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
        rel_abunds = abundance_trajectories[:, tidx] / total_mass[tidx]  # indexed by taxa

        phis = rel_abunds * total_reads
        epsilons = negbin_a1 + (negbin_a0 / rel_abunds)
        reads_all.append(
            np.array([
                negative_binomial(phi=phi, eps=eps, rng=rng)
                if ~np.isinf(eps)
                else 0
                for phi, eps in zip(phis, epsilons)
            ])
        )
    return np.stack(reads_all, axis=1)


def negative_binomial(phi, eps, rng: np.random.Generator) -> int:
    """
    :param phi: the mean parameter
    :param eps: the dispersion parameter
    :return:

    Refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html
    which gives:
     p = (mean) / (var) = 1 / (1 + phi * eps)
     n = (mean^2) / (var - mean) = 1 / eps
    """
    p = 1 / (1 + phi * eps)
    n = 1 / eps
    return rng.negative_binomial(n=n, p=p)


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


# ====== Dataframe formatters
def name_formatter(subj_idx: int, t_idx: int) -> str:
    return f"{subj_idx}-T{t_idx}"


def create_metadata_df(n_subjects: int, timepoints: List[float]) -> pd.DataFrame:
    df_entries = []
    for subj_idx in range(n_subjects):
        subj_name = generate_mouse_name(subj_idx)
        for t_idx, t in enumerate(timepoints):
            sample_id = name_formatter(subj_idx, t_idx)
            df_entries.append({
                "sampleID": sample_id,
                "subject": subj_name,
                "time": t
            })
    return pd.DataFrame(df_entries)


def create_counts_df(
        n_subjects: int,
        timepoints: List[float],
        read_counts: Dict[str, np.ndarray],
        taxa: TaxaSet,
) -> pd.DataFrame:
    df_entries = [
        {
            "taxaName": asv.name,
            "sampleID": name_formatter(subj_idx, t_idx),
            "reads": read_counts[generate_mouse_name(subj_idx)][asv_index, t_idx]
        }
        for asv_index, asv in enumerate(taxa)
        for subj_idx in range(n_subjects)
        for t_idx, _ in enumerate(timepoints)
    ]
    return pd.DataFrame(df_entries)


def create_qpcr_df(
        n_subjects: int,
        timepoints: List[float],
        qpcr_values: Dict[str, np.ndarray],
) -> pd.DataFrame:
    df_entries = [
        {
            "sampleID": name_formatter(subj_idx, t_idx),
            "measurement1": qpcr_values[generate_mouse_name(subj_idx)][t_idx, 0],
            "measurement2": qpcr_values[generate_mouse_name(subj_idx)][t_idx, 1],
            "measurement3": qpcr_values[generate_mouse_name(subj_idx)][t_idx, 2]
        }
        for subj_idx in range(n_subjects)
        for t_idx, _ in enumerate(timepoints)
    ]
    return pd.DataFrame(df_entries)


def create_replicate_metadata_df(n_replicates: int, integer_timepoints: List[int]) -> pd.DataFrame:
    df_entries = []
    for t in integer_timepoints:
        for replicate_idx in range(n_replicates):
            df_entries.append({
                "sampleID": f'M0-D{t}-R{replicate_idx}',
                "subject": f'M0-D{t}',
                "time": replicate_idx  # this looks wrong at first glance, but the point is that each physical replicate is treated as a separate sample.
            })
    return pd.DataFrame(df_entries)


def create_replicate_counts_df(
        n_replicates: int,
        integer_timepoints: List[int],
        read_counts: np.ndarray,
        taxa: TaxaSet,
) -> pd.DataFrame:
    df_entries = [
        {
            "taxaName": asv.name,
            "sampleID": f'M0-D{t}-R{replicate_idx}',
            "reads": read_counts[replicate_idx][asv_index, t_idx]
        }
        for asv_index, asv in enumerate(taxa)
        for replicate_idx in range(n_replicates)
        for t_idx, t in enumerate(integer_timepoints)
    ]
    return pd.DataFrame(df_entries)


def create_replicate_qpcr_df(
        n_replicates: int,
        integer_timepoints: List[int],
        qpcr_values: np.ndarray,
) -> pd.DataFrame:
    df_entries = [
        {
            "sampleID": f'M0-D{t}-R{replicate_idx}',
            "measurement1": qpcr_values[replicate_idx][t_idx, 0],
            "measurement2": qpcr_values[replicate_idx][t_idx, 1],
            "measurement3": qpcr_values[replicate_idx][t_idx, 2]
        }
        for replicate_idx in range(n_replicates)
        for t_idx, t in enumerate(integer_timepoints)
    ]
    return pd.DataFrame(df_entries)


# ====== Main sampling functions
def create_synthetic_dataset(
        source_taxa: md2.TaxaSet,
        forward_sims: np.ndarray,
        sim_timepoints: np.ndarray,
        sim_read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    read_samples_per_subj, qpcr_samples_per_subj, subset_timepoints = sample_data_from_fwsim(
        forward_sims,
        [generate_mouse_name(i) for i in range(forward_sims.shape[0])],
        sim_timepoints,
        sim_read_depth,
        negbin_a0,
        negbin_a1,
        qpcr_noise_scale,
        rng,
    )

    n_subjects = len(read_samples_per_subj)
    metadata_df = create_metadata_df(n_subjects, subset_timepoints)
    counts_df = create_counts_df(n_subjects, subset_timepoints, read_samples_per_subj, source_taxa)
    qpcr_df = create_qpcr_df(n_subjects, subset_timepoints, qpcr_samples_per_subj)
    return counts_df, qpcr_df, metadata_df


def create_synthetic_replicates(
        source_taxa: md2.TaxaSet,
        forward_sims: np.ndarray,
        sim_timepoints: np.ndarray,
        source_subj_idx: int,
        negbin_a0: float,
        negbin_a1: float,
        num_reads: int,
        qpcr_noise_scale: float,
        n_physical_replicates: int,
        rng: np.random.Generator
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Extract the "true" biomass (the posterior median of the filtered state abundances)
    read_samples_per_replicate, qpcr_samples_per_replicate, integer_timepoints = sample_replicate_data_from_fwsim(
        forward_sims[source_subj_idx],
        sim_timepoints,
        n_physical_replicates,
        num_reads,
        negbin_a0,
        negbin_a1,
        qpcr_noise_scale,
        rng,
    )

    metadata_df = create_replicate_metadata_df(n_physical_replicates, integer_timepoints)
    counts_df = create_replicate_counts_df(n_physical_replicates, integer_timepoints, read_samples_per_replicate, source_taxa)
    qpcr_df = create_replicate_qpcr_df(n_physical_replicates, integer_timepoints, qpcr_samples_per_replicate)
    return counts_df, qpcr_df, metadata_df


def main(
        source_study: md2.Study,
        trajectories: np.ndarray,
        timepoints: np.ndarray,
        out_dir: Path,
        sim_read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        calibration_physical_replicates: int,
        rng: np.random.Generator
):
    out_dir.mkdir(exist_ok=True, parents=True)
    assert len(timepoints) == trajectories.shape[-1]

    # n_subj, n_taxa, n_timepoints = trajectories.shape
    inference_counts_df, inference_qpcr_df, inference_metadata_df = create_synthetic_dataset(
        source_study.taxa,
        trajectories, timepoints,
        sim_read_depth,
        negbin_a0, negbin_a1,
        qpcr_noise_scale,
        rng,
    )
    replicate_counts_df, replicate_qpcr_df, replicate_metadata_df = create_synthetic_replicates(
        source_study.taxa,
        trajectories, timepoints,
        0,
        negbin_a0, negbin_a1,
        sim_read_depth,
        qpcr_noise_scale,
        calibration_physical_replicates,
        rng,
    )
    inference_counts_df.to_csv(out_dir / "counts.tsv", sep='\t', index=False)
    inference_qpcr_df.to_csv(out_dir / "qpcr.tsv", sep='\t', index=False)
    inference_metadata_df.to_csv(out_dir/ "metadata.tsv", sep='\t', index=False)
    replicate_counts_df.to_csv(out_dir / "replicate_counts.tsv", sep='\t', index=False)
    replicate_qpcr_df.to_csv(out_dir / "replicate_qpcr.tsv", sep='\t', index=False)
    replicate_metadata_df.to_csv(out_dir / "replicate_metadata.tsv", sep='\t', index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-study", "-i", dest="source_study_path", required=True)
    parser.add_argument("--sim-dir", "-sd", dest="sim_dir", type=str, required=True),
    parser.add_argument("--outdir", "-o", dest="out_dir", type=str, required=True),
    parser.add_argument("--read-depth", "-r", dest="read_depth", type=int, required=True)
    parser.add_argument("--a0", "-a0", dest="negbin_a0", type=float, required=True)
    parser.add_argument("--a1", "-a1", dest="negbin_a1", type=float, required=True)
    parser.add_argument("--qpcr-noise-scale", "-q", dest="qpcr_noise_scale", type=float, required=True)
    parser.add_argument(
        "--calibration-replicates", "-c", dest="calibration_physical_replicates",
        type=int, required=True,
        help='The number of physical fecal replicates to simulate per sample. Used for the replicate dataset.'
    )
    parser.add_argument("--seed", "-s", dest="seed", type=int, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    sim_dir = Path(args.sim_dir)
    main(
        source_study=md2.Study.load(args.source_study_path),
        trajectories=np.load(sim_dir / "trajectories.npy"),
        timepoints=np.load(sim_dir / "timepoints.npy"),
        out_dir=Path(args.out_dir),
        sim_read_depth=args.read_depth,
        negbin_a0=args.negbin_a0,
        negbin_a1=args.negbin_a1,
        qpcr_noise_scale=args.qpcr_noise_scale,
        calibration_physical_replicates=args.calibration_physical_replicates,
        rng=np.random.default_rng(args.seed)
    )