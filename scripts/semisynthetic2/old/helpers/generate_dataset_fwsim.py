from pathlib import Path
from typing import Tuple, Dict, List
import pickle as pkl

import argparse
import numpy as np

import mdsine2 as md2
from mdsine2 import TaxaSet, Perturbations, BasePerturbation
from scripts.semisynthetic2.data_generation.python_helpers.base import forward_simulate_subject


# ======= Helper functions
# the timepoints, perturbation windows are the same for all simulated subjects.
def create_synthetic_dataset(
        source_taxa: TaxaSet,
        n_subjects: int,
        forward_sims: Dict[str, np.ndarray],
        timepoints: List[float],
        pert_names: List[str],
        pert_starts: List[float],
        pert_ends: List[float],
        sim_read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator,
) -> md2.Study:
    print("Sampling.")
    assert len(pert_names) == len(pert_starts)
    assert len(pert_names) == len(pert_ends)

    read_samples_per_subj, qpcr_samples_per_subj = sample_data_from_fwsim(
        source_taxa,
        n_subjects,
        timepoints,
        forward_sims,
        sim_read_depth,
        negbin_a0,
        negbin_a1,
        qpcr_noise_scale,
        rng,
    )

    # create the study object.
    synthetic_study = md2.Study(taxa=source_taxa, name='simulated')
    for subj_idx in range(n_subjects):
        subj_name = generate_subject_name(subj_idx)
        synthetic_study.add_subject(name=subj_name)
        new_subj = synthetic_study[subj_name]

        for tidx, t in enumerate(timepoints):
            sampled_reads = read_samples_per_subj[subj_name]  # (n_taxa, n_timepoints)
            sampled_qpcr = qpcr_samples_per_subj[subj_name]
            new_subj.reads[t] = sampled_reads[:, tidx]
            new_subj.qpcr[t] = md2.qPCRdata(cfus=sampled_qpcr[tidx], mass=1., dilution_factor=1.)
        new_subj.times = timepoints

    perts_obj = Perturbations()
    synthetic_study.perturbations = perts_obj
    for pert_name, pert_start_t, pert_end_t in zip(pert_names, pert_starts, pert_ends):
        perts_obj.append(BasePerturbation(
            name=pert_name,
            starts={generate_subject_name(subj_idx): pert_start_t for subj_idx in range(n_subjects)},
            ends={generate_subject_name(subj_idx): pert_end_t for subj_idx in range(n_subjects)}
        ))
    print("synthetic study name: {}".format(synthetic_study.name))
    return synthetic_study


def create_synthetic_replicates(
        source_study: md2.Study,
        forward_sims: Dict[str, np.ndarray],
        source_subj_name: str,
        source_subj_timepoints: List[float],  # timepoints
        num_physical_replicates: int,
        replicate_study_name: str,
        negbin_a0: float,
        negbin_a1: float,
        num_reads: int,
        qpcr_noise_scale: float,
        rng: np.random.Generator
) -> md2.Study:
    # Extract the source subj and related data
    src_subj = source_study[source_subj_name]
    timepoint_indices = {t: t_idx for t_idx, t in enumerate(src_subj.times)}

    # Extract the "true" biomass (the posterior median of the filtered state abundances)
    true_traj = forward_sims[source_subj_name]

    # Make the study object
    replicate_study = md2.Study(taxa=source_study.taxa, name=replicate_study_name)
    for t in source_subj_timepoints:
        replicate_study.add_subject(name=f'M{source_subj_name}-D{t}')

    # Add times for each subject
    replicate_times = np.arange(0., num_physical_replicates, 1.0)
    for replicate_subj in replicate_study:
        replicate_subj.times = replicate_times

    for replicate_subj, subj_t in zip(replicate_study, source_subj_timepoints):
        # Extract timepoint abundances
        if subj_t not in timepoint_indices:
            raise Exception("Couldn't find timepoint {} in subject {}.".format(subj_t, source_subj_name))
        t_idx = timepoint_indices[subj_t]

        conc = true_traj[:, t_idx]
        total_mass = np.sum(conc)
        rel_mass = conc / total_mass

        # simulate physical replicates
        for repl_t in replicate_times:
            # Make the reads
            replicate_subj.reads[repl_t] = np.array([
                negative_binomial(phi=ra_ratio * num_reads, eps=negbin_a1 + (negbin_a0 / ra_ratio), rng=rng)
                if ra_ratio > 0
                else 0
                for ra_ratio in rel_mass
            ])

            # simulate technical/measurement triplicates
            triplicates = np.exp(
                np.log(total_mass) +
                rng.normal(loc=0.0, scale=qpcr_noise_scale, size=3)
            )
            replicate_subj.qpcr[repl_t] = md2.qPCRdata(cfus=triplicates, mass=1., dilution_factor=1.)
    print("replicate study name: {}".format(replicate_study.name))
    return replicate_study


def generate_subject_name(subj_idx: int) -> str:
    return  "SUBJ_{}".format(subj_idx)


def sample_data_from_fwsim(
        source_taxa: TaxaSet,
        n_subjects: int,
        timepoints: List[float],
        forward_sims: Dict[str, np.ndarray],
        read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    :return: A tuple.
    (1) A dictionary of (Taxa x Timepoint) abundance samples per subject,
    (2) A dictionary of (Timepoint x Replicate) qPCR samples per subject.
    """
    reads_per_subject = {}
    qpcr_per_subject = {}
    for subj_idx in range(n_subjects):
        subj_name = generate_subject_name(subj_idx)
        trajs = forward_sims[subj_name]

        if trajs.shape[0] != len(source_taxa):
            raise ValueError("Matrix's n_taxa dimension ({}) doesn't match taxaset size ({})".format(
                trajs.shape[0], len(source_taxa)
            ))
        if trajs.shape[1] != len(timepoints):
            raise ValueError("Matrix's n_times dimension ({}) doesn't match subj times size ({})".format(
                trajs.shape[1], len(timepoints)
            ))

        read_counts = sample_reads(read_depth, trajs, negbin_a0, negbin_a1, rng=rng)  # n_taxa x n_timepoints
        reads_per_subject[subj_name] = read_counts

        qpcrs = sample_qpcr(np.sum(trajs, axis=0), qpcr_noise_scale, n_replicates=3, rng=rng)
        qpcr_per_subject[subj_name] = qpcrs
    return reads_per_subject, qpcr_per_subject


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


# ============ CLI interface
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--study', '-s', dest="study_path",
        type=str, required=True,
        help="The path to the original (real) datas Study pickle file."
    )
    parser.add_argument(
        '--truth-dir', '-t', dest='ground_truth_dir',
        type=str, required=True,
        help="A directory to store ground truth values."
    )
    parser.add_argument(
        '--out', '-o', dest='out_path',
        type=str, required=True,
        help="The desired full path to which the synthetic study object will be saved."
    )
    parser.add_argument(
        '--replicate-out', '-ro', dest='replicate_out_path',
        type=str, required=True,
        help="The desired full path to which the synthetic study object will be saved.."
    )
    parser.add_argument(
        '--num-subjects', '-ns', dest='num_subjects',
        type=int, required=True,
        help="The number of subjects to simulate."
    )
    parser.add_argument(
        '--read-depth', '-r', dest='read_depth',
        type=int, required=True,
        help="The overall read depth to simulate per timepoint."
    )
    parser.add_argument(
        '--a0', '-a0', dest='negbin_a0',
        type=float, required=True,
        help="The a0 parameter of the negbin parametrization."
    )
    parser.add_argument(
        '--a1', '-a1', dest='negbin_a1',
        type=float, required=True,
        help="The a1 parameter of the negbin parametrization."
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
    parser.add_argument('--sim-dt', type=float, required=False, default=0.01)
    parser.add_argument('--sim-max', type=float, required=False, default=1e20)
    return parser.parse_args()


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


def main(
        study_path: Path,
        ground_truth_dir: Path,
        out_path: Path,
        num_subjects: int,
        replicate_out_path: Path,
        read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        seed: int,
        sim_max: float,
        sim_dt: float
):
    source_study = md2.Study.load(str(study_path))
    ground_truth_dir.mkdir(parents=True, exist_ok=True)
    print(f"Using ground truth dir {ground_truth_dir}")

    # Load ground truth parameters.
    glv_sim_path = ground_truth_dir / 'glv_best_sim.pkl'
    with open(glv_sim_path, 'rb') as f:
        glv_params = pkl.load(f)

    # Initialize.
    rng = np.random.default_rng(seed)
    out_path.parent.mkdir(exist_ok=True, parents=True)

    # Get some metadata (perturbation names and times)
    target_subj = subj_with_most_timepoints(source_study)  # pick subject with most # of timepoints

    pert_names = [p.name for p in source_study.perturbations]
    pert_starts = [
        pert.starts[target_subj.name]
        for pert in source_study.perturbations
    ]
    pert_ends = [
        pert.ends[target_subj.name]
        for pert in source_study.perturbations
    ]

    timepoints = target_subj.times

    n_taxa = len(source_study.taxa)
    initial_conditions = np.random.randn(size=(num_subjects, n_taxa))

    # Simulate all the ground-truth trajectories.
    forward_sims = {}
    for subj_idx in range(num_subjects):
        sims, times = forward_simulate_subject(
            glv_params,
            pert_starts,
            pert_ends,
            timepoints,
            initial_conditions[subj_idx],
            sim_dt,
            sim_max,
            time_subset=True
        )
        assert len(times) == len(timepoints)
        assert np.sum(times != timepoints) == 0
        forward_sims[generate_subject_name(subj_idx)] = sims

    synth_study = create_synthetic_dataset(
        source_taxa=source_study.taxa,
        n_subjects=num_subjects,
        forward_sims=forward_sims,
        timepoints=timepoints,
        pert_names=pert_names,
        pert_starts=pert_starts,
        pert_ends=pert_ends,
        sim_read_depth=read_depth,
        negbin_a0=negbin_a0,
        negbin_a1=negbin_a1,
        qpcr_noise_scale=qpcr_noise_scale,
        rng=rng,
    )

    synth_study.save(str(out_path))
    print(f"Saved semisynthetic dataset to {out_path}.")

    replicate_timepoints = timepoints[6:6+3]  # somewhat arbitrary, but similar to real data
    synth_repl_study = create_synthetic_replicates(
        source_study=synth_study,
        forward_sims=forward_sims,
        source_subj_name=next(iter(synth_study)).name,
        source_subj_timepoints=replicate_timepoints,  # timepoints from which to source the samples
        num_physical_replicates=6,
        replicate_study_name='synthetic-replicate',
        negbin_a0=negbin_a0,
        negbin_a1=negbin_a1,
        num_reads=read_depth,
        qpcr_noise_scale=qpcr_noise_scale,
        rng=rng
    )
    replicate_out_path.parent.mkdir(exist_ok=True, parents=True)
    synth_repl_study.save(str(replicate_out_path))
    print(f"Saved replicates to {replicate_out_path}.")


if __name__ == "__main__":
    args = parse_args()
    main(
        study_path=Path(args.study_path),
        ground_truth_dir=Path(args.ground_truth_dir),
        num_subjects=args.num_subjects,
        out_path=Path(args.out_path),
        replicate_out_path=Path(args.replicate_out_path),
        read_depth=args.read_depth,
        negbin_a0=args.negbin_a0,
        negbin_a1=args.negbin_a1,
        qpcr_noise_scale=args.qpcr_noise_scale,
        seed=args.seed,
        sim_max=args.sim_max,
        sim_dt=args.sim_dt
    )
