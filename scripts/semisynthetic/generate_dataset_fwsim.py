from collections import namedtuple
from pathlib import Path
from typing import Tuple, Dict, List, Set
import pickle as pkl

import argparse
import numpy as np
from tqdm import tqdm

import mdsine2 as md2
import matplotlib.pyplot as plt

from mdsine2 import TaxaSet


# ======= Helper functions
def create_synthetic_dataset(
        source_study: md2.Study,
        forward_sims: Dict[str, np.ndarray],
        sim_read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator,
) -> md2.Study:
    print("Sampling.")
    read_samples_per_subj, qpcr_samples_per_subj = sample_data_from_fwsim(
        source_study,
        forward_sims,
        sim_read_depth,
        negbin_a0,
        negbin_a1,
        qpcr_noise_scale,
        rng,
    )

    # # Prune all taxa that don't have some nonzero read. Ensure that the taxa slice is deleted from reads and qpcr also.
    # print("Pruning.")
    # concat_reads_count = np.hstack([subj_reads for subj_name, subj_reads in read_samples_per_subj.items()]).sum(axis=-1)
    # taxa_idx_to_keep, = np.where(concat_reads_count > 0)
    # taxa_idx_set = set(taxa_idx_to_keep)
    #
    # new_taxa = TaxaSet()  # create new taxaset
    # _i = 0
    # for taxa_idx, taxon in enumerate(source_study.taxa):
    #     if taxa_idx in taxa_idx_set:
    #         new_taxa.add(taxon)
    #         taxon.idx = _i
    #         _i += 1
    #
    # read_samples_per_subj = {
    #     subj_name: subj_reads[taxa_idx_to_keep, :]  # prune
    #     for subj_name, subj_reads in read_samples_per_subj.items()
    # }
    #
    # print("Retained {} of {} taxa.".format(
    #     taxa_idx_to_keep.shape[0],
    #     len(source_study.taxa)
    # ))

    # create the study object.
    synthetic_study = md2.Study(taxa=source_study.taxa)
    for source_subj in source_study:
        synthetic_study.add_subject(name=source_subj.name)
        new_subj = synthetic_study[source_subj.name]
        for tidx, t in enumerate(source_subj.times):
            sampled_reads = read_samples_per_subj[source_subj.name]  # (n_taxa, n_timepoints)
            sampled_qpcr = qpcr_samples_per_subj[source_subj.name]
            new_subj.reads[t] = sampled_reads[:, tidx]
            new_subj.qpcr[t] = md2.qPCRdata(cfus=sampled_qpcr[tidx], mass=1., dilution_factor=1.)
        new_subj.times = source_subj.times
    synthetic_study.perturbations = source_study.perturbations
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
    return replicate_study



GLVParamSet = namedtuple(
    'GLVParamSet',
    [
        'growth', 'interactions', 'perturbations'
    ]
)


def extract_glv_model(
        study: md2.Study,
        growths: np.ndarray,
        interactions: np.ndarray,
        perturbations: List[np.ndarray],
        coclusterings: np.ndarray,
        dt: float,
        sim_max: float,
) -> Tuple[GLVParamSet, np.ndarray, Dict[str, np.ndarray]]:
    """Pick a best forward simulation based on error evaluation metric."""

    subsample_every = 1000
    if subsample_every > 1:
        print("Subsampling every {} MCMC samples. (overall, there are {} samples)".format(subsample_every, interactions.shape[0]))
    params = [
        GLVParamSet(
            growth=growths[i],
            interactions=interactions[i],
            perturbations=[p[i] for p in perturbations]
        )
        for i in range(0, interactions.shape[0], subsample_every)
    ]
    coclusterings = coclusterings[::subsample_every]

    # Extract the best parameter set using forward simulations.
    param_fwsim_errors = []
    param_num_surviving_modules = []
    # for param_set, coclust_mat in tqdm(zip(params, coclusterings), total=len(params)):
    iter_idx = 0
    for param_set, coclust_mat in zip(params, coclusterings):
        err, surviving_taxa, surviving_modules, total_modules = evaluate_parameter_fwsim(param_set, study, dt, sim_max, coclust_mat)
        if len(surviving_taxa) < len(study.taxa):
            param_fwsim_errors.append(np.inf)
        else:
            param_fwsim_errors.append(err)
        param_num_surviving_modules.append(len(surviving_modules))
        print("[iter={}] Modules survived: {} / {}".format(iter_idx, len(surviving_modules), total_modules))
        iter_idx += subsample_every

    best_idx = np.argmin(param_fwsim_errors)
    if param_fwsim_errors[best_idx] == np.inf:
        raise ValueError("Couldn't find best sim with all surviving taxa")
    print("Choice: index {} (total {} entries)".format(
        best_idx,
        len(params)
    ))
    return params[best_idx], coclusterings[best_idx], forward_simulate(params[best_idx], study, dt, sim_max)


def forward_simulate(
        glv_params: GLVParamSet,
        study: md2.Study,
        dt: float,
        sim_max: float
) -> Dict[str, np.ndarray]:
    """Forward simulation for all subjects in a Study object"""
    return {
        subj.name: forward_simulate_subject(glv_params, study, subj, dt, sim_max)[0]
        for subj in study
    }


def forward_simulate_subject(
        glv_params: GLVParamSet,
        study: md2.Study,
        subject: md2.Subject,
        dt: float,
        sim_max: float,
        time_subset: bool = True
) -> np.ndarray:
    """Forward simulation for a single subject"""
    # ======= Perturbations
    perturbations_start = []
    perturbations_end = []
    if study.perturbations is not None:
        for pert in study.perturbations:
            perturbations_start.append(pert.starts[subject.name])
            perturbations_end.append(pert.ends[subject.name])

    # print("Growth:", glv_params.growth)
    # print("Diag:", np.diag(glv_params.interactions))
    # print("Start day: {}".format(subject.times[0]))
    dyn = md2.model.gLVDynamicsSingleClustering(
        growth=glv_params.growth,
        interactions=glv_params.interactions,
        perturbations=glv_params.perturbations,
        perturbation_starts=perturbations_start,
        perturbation_ends=perturbations_end,
        start_day=subject.times[0],
        sim_max=sim_max
    )

    # print("Data shape: ", subject.matrix()['abs'].shape)
    initial_conditions = subject.matrix()['abs'][:, 0] + 1e4
    if time_subset:
        x = md2.integrate(
            dynamics=dyn,
            initial_conditions=np.expand_dims(initial_conditions, 1),
            dt=dt,
            final_day=subject.times[-1],
            subsample=True,
            times=subject.times
        )
    else:
        x = md2.integrate(
            dynamics=dyn,
            initial_conditions=np.expand_dims(initial_conditions, 1),
            dt=dt,
            final_day=subject.times[-1],
            subsample=False
        )
    fwsim_values = x['X']
    times = x['times']
    return fwsim_values, times


# =========== helpers
def extract_clustering_from_matrix(_mat) -> np.ndarray:
    n_items = _mat.shape[0]
    items_left = set(range(n_items))

    clustering_assignments = np.zeros(n_items, dtype=int)

    c_idx = -1
    while len(items_left) > 0:
        c_idx += 1  # new cluster
        x = items_left.pop()  # member
        clustering_assignments[x] = c_idx

        # get all guys in the same cluster
        to_remove = set()
        for y in items_left:
            if _mat[x,y]:
                clustering_assignments[y] = c_idx
                to_remove.add(y)
        items_left = items_left.difference(to_remove)
    return clustering_assignments


def evaluate_parameter_fwsim(
        params: GLVParamSet,
        study: md2.Study,
        dt: float,
        sim_max: float,
        coclust_matrix: np.ndarray
) -> Tuple[float, Set[int], Set[int], int]:
    """
    Evaluate the fwsim error for each subject (And sum them).
    """
    fwsims = forward_simulate(params, study, dt, sim_max)
    errors = []
    eps = 1e-5
    for subj_name, fwsim_trajs in fwsims.items():
        subj = study[subj_name]
        measurements = subj.matrix()['abs']

        # RMS-log10 error, excluding the points where data says zero.
        # fwsim_trajs: [taxa x timepoints] array
        # measurements: [taxa x timepoints] array, with zeroes in some of the entries

        mask = measurements > 0
        diff_logs = np.log10(fwsim_trajs + eps) - np.log10(measurements + eps)
        subj_err = np.sqrt(np.mean(np.square(diff_logs[mask])))
        errors.append(subj_err)

    # What taxa survived? (Only keep taxa which surpasses 1e5 abundance for some timepoint in some synthetic mouse)
    taxa_max_abundance = np.stack([
        fwsim_trajs.max(axis=-1)  # length=n_taxa, after maximizing across timepoints
        for subj_name, fwsim_trajs in fwsims.items()
    ], axis=0).max(axis=0)  # length=n_taxa, max across subjects

    clust_assignments = extract_clustering_from_matrix(coclust_matrix)
    leftover_module_set = set()
    leftover_taxa_set = set()
    for taxa_idx, taxa in enumerate(study.taxa):
        if taxa_max_abundance[taxa_idx] > 1e5:
            leftover_module_set.add(clust_assignments[taxa_idx])
            leftover_taxa_set.add(taxa_idx)
    total_modules = np.max(clust_assignments) + 1

    print("leftover taxa: {}".format(len(leftover_taxa_set)))
    return np.sum(errors), leftover_taxa_set, leftover_module_set, total_modules


def sample_data_from_fwsim(
        study: md2.Study,
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
    for subj_idx, subj in enumerate(study):
        trajs = forward_sims[subj.name]

        if trajs.shape[0] != len(study.taxa):
            raise ValueError("Matrix's n_taxa dimension ({}) doesn't match taxaset size ({})".format(
                trajs.shape[0], len(study.taxa)
            ))
        if trajs.shape[1] != len(subj.times):
            raise ValueError("Matrix's n_times dimension ({}) doesn't match subj times size ({})".format(
                trajs.shape[1], len(subj.times)
            ))

        read_counts = sample_reads(read_depth, trajs, negbin_a0, negbin_a1, rng=rng)  # n_taxa x n_timepoints
        reads_per_subject[subj.name] = read_counts

        qpcrs = sample_qpcr(np.sum(trajs, axis=0), qpcr_noise_scale, n_replicates=3, rng=rng)
        qpcr_per_subject[subj.name] = qpcrs
    return reads_per_subject, qpcr_per_subject


def render_plots(study: md2.Study, sims: Dict[str, np.ndarray], plot_dir: Path):
    width = 6
    height_per_subj = 6
    fmt = 'pdf'
    for taxa_idx, taxa in enumerate(study.taxa):
        plot_path = plot_dir / f'{taxa.name}.{fmt}'

        fig, axes = plt.subplots(len(study), 1, figsize=(width, height_per_subj * len(study)))
        for subj, ax in zip(study, axes):
            measurements = subj.matrix()['abs'][taxa_idx, :]
            trajs = sims[subj.name]
            ax.plot(subj.times, trajs[taxa_idx, :])
            ax.plot(subj.times, measurements, marker='x', color='black', linestyle=':')
            ax.set_title(f'Subject {subj.name}')
            ax.set_yscale('log')
        fig.suptitle(f"Taxa {taxa.name}")
        plt.savefig(plot_path, format=fmt)
        plt.close(fig)


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
        '--growths', dest='growth_path',
        type=str, required=True, help="The path to the .npy growth rate MCMC samples."
    )
    parser.add_argument(
        '--interactions', dest='interactions_path',
        type=str, required=True, help="The path to the .npy interactions MCMC samples."
    )
    parser.add_argument(
        '--perts', dest='perts_path',
        type=str, required=True, help="The path to the .npz perturbation strengths MCMC samples."
    )
    parser.add_argument(
        '--coclust', dest='coclust_path',
        type=str, required=True, help="The path to the .npz coclustering MCMC samples."
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


def main(
        study_path: Path,
        growth_path: Path,
        interactions_path: Path,
        perturbations_path: Path,
        coclusterings_path: Path,
        ground_truth_dir: Path,
        out_path: Path,
        replicate_out_path: Path,
        read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        seed: int,
        sim_max: float,
        sim_dt: float
):
    """
    Read the fixed-module inference, use the stored "filtered-state" (latent traj X) samples.

    :param fixed_module_pkl_path:
    :param read_depth:
    :return:
    """
    source_study = md2.Study.load(str(study_path))
    ground_truth_dir.mkdir(parents=True, exist_ok=True)
    print(f"Using ground truth dir {ground_truth_dir}")

    glv_sim_path = ground_truth_dir / 'glv_best_sim.pkl'
    if glv_sim_path.exists():
        with open(glv_sim_path, 'rb') as f:
            glv_params, forward_sims = pkl.load(f)
    else:
        # Take the median parameter set from posterior.
        growths = np.load(growth_path)
        interactions = np.load(interactions_path)
        perturbations_map = np.load(perturbations_path)
        perturbations = [
            perturbations_map[pert.name]
            for pert in source_study.perturbations
        ]
        coclusterings = np.load(coclusterings_path)
        assert coclusterings.shape[0] == interactions.shape[0]  # ensure that the number of samples match

        glv_params, instance_coclusters, forward_sims = extract_glv_model(source_study, growths, interactions, perturbations, coclusterings, sim_dt, sim_max)
        plot_dir = ground_truth_dir / 'fwsim-plots'
        plot_dir.mkdir(exist_ok=True, parents=True)
        render_plots(source_study, forward_sims, plot_dir)
        with open(glv_sim_path, 'wb') as f:
            pkl.dump((glv_params, forward_sims), f)
        np.save(ground_truth_dir / "coclusters.npy", instance_coclusters)

        # Save the full forward simulation.
        for subj in source_study:
            sims, times = forward_simulate_subject(
                glv_params,
                source_study,
                subj,
                sim_dt,
                sim_max,
                time_subset=False
            )
            np.savez(ground_truth_dir / f'forward_sim_full_{subj.name}.npz', sims=sims, times=times)

    rng = np.random.default_rng(seed)
    out_path.parent.mkdir(exist_ok=True, parents=True)
    synth_study = create_synthetic_dataset(
        source_study=source_study,
        forward_sims=forward_sims,
        sim_read_depth=read_depth,
        negbin_a0=negbin_a0,
        negbin_a1=negbin_a1,
        qpcr_noise_scale=qpcr_noise_scale,
        rng=rng,
    )

    synth_study.save(str(out_path))
    print(f"Saved semisynthetic dataset to {out_path}.")

    synth_repl_study = create_synthetic_replicates(
        source_study=source_study,
        forward_sims=forward_sims,
        source_subj_name='2',
        source_subj_timepoints=[8.0, 9.0, 10.0],  # timepoints
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
        growth_path=Path(args.growth_path),
        interactions_path=Path(args.interactions_path),
        perturbations_path=Path(args.perts_path),
        coclusterings_path=Path(args.coclust_path),
        ground_truth_dir=Path(args.ground_truth_dir),
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
