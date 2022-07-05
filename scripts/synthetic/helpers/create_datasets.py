"""
Python script for generating semisynthetic samples for a given seed + noise level.
Takes as input MDSINE1's BVS sample matrix file
"""
from typing import Tuple, List
from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np

from mdsine2 import *
from mdsine2.names import STRNAMES


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_glv_params', type=str, required=True,
                        help='<Required> The path to a .npz file containing `growth_rates`, `interactions` arrays')
    parser.add_argument('-t', '--time_points_file', type=str, required=True,
                        help='<Required> A path to a text file containing a list of time points to pull out gLV '
                             'measurements from.')
    parser.add_argument('-n', '--num_subjects', type=int, required=True,
                        help='<Required> The number of subjecs to simulate to lump into a single cohort.')
    parser.add_argument('-o', '--out_dir', type=str, required=True,
                        help='<Required> The directory to output the sampled subject to.')
    parser.add_argument('-s', '--seed', type=int, required=True,
                        help='<Required> The seed to use for random sampling.')

    # Optional parameters
    parser.add_argument('-p', '--process_var', type=float, required=False, default=0.01)
    parser.add_argument('-dt', '--sim_dt', type=float, required=False, default=0.01)
    parser.add_argument('-a', '--dmd_alpha_scale', type=float, required=False, default=286,
                        help='DMD dispersion parameter estimated from mean estimates at all time-points from '
                             'C. diff data (Default: 286, carry-over from MDSINE1)')
    parser.add_argument('-r', '--read_depth', type=int, required=False, default=50000)
    parser.add_argument('--low_noise', type=float, required=False, default=0.01)
    parser.add_argument('--medium_noise', type=float, required=False, default=0.1)
    parser.add_argument('--high_noise', type=float, required=False, default=0.2)
    return parser.parse_args()


def make_synthetic(
        name: str,
        taxa: TaxaSet,
        growth_rate_values: np.ndarray,
        interaction_values: np.ndarray,
        interaction_indicators: np.ndarray,
        seed: int
) -> Synthetic:
    syn = Synthetic(name=name, seed=seed)
    syn.taxa = taxa

    clustering = Clustering(clusters=None, G=syn.G, items=syn.taxa, name=STRNAMES.CLUSTERING_OBJ)
    interactions = Interactions(clustering=clustering, use_indicators=True, name=STRNAMES.INTERACTIONS_OBJ, G=syn.G)
    for interaction in interactions:
        # Set interaction values
        target_cid = interaction.target_cid
        source_cid = interaction.source_cid

        tcidx = clustering.cid2cidx[target_cid]
        scidx = clustering.cid2cidx[source_cid]

        interaction.value = interaction_values[tcidx, scidx]
        interaction.indicator = interaction_indicators[tcidx, scidx]

    syn.model.interactions = interaction_values
    syn.model.growth = growth_rate_values
    return syn


def parse_glv_params(params_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], np.ndarray, np.ndarray]:
    params = np.load(str(params_path))
    growth = params['growth_rates']
    interactions = params['interactions']
    initial_cond_mean = params['initial_mean'] * np.ones(len(growth), dtype=float)
    initial_cond_std = params['initial_std'] * np.ones(len(growth), dtype=float)
    indicators = (interactions != 0.0)
    taxa_names = [f'TAXA_{i+1}' for i in range(len(growth))]
    return growth, interactions, indicators, taxa_names, initial_cond_mean, initial_cond_std


def parse_time_points(time_points_path: Path) -> np.ndarray:
    time_points = []
    with open(time_points_path, 'rt') as f:
        for line in f:
            t = float(line.rstrip())
            time_points.append(t)
    return np.array(time_points, dtype=float)


def simulate_reads_dmd(synth: Synthetic, study_name: str, alpha_scale: float, num_reads: int, qpcr_noise_scale: float) -> Study:
    # Make the study object
    study = Study(taxa=synth.taxa, name=study_name)
    for subjname in synth.subjs:
        study.add_subject(name=subjname)

    # Add times for each subject
    for subj in study:
        subj.times = synth.times

    for subj in study:
        total_mass = np.sum(synth._data[subj.name], axis=0)  # length T
        for tidx, t in enumerate(subj.times):
            # Make the reads
            rel_abund = synth._data[subj.name][:, tidx] / total_mass[tidx]
            alpha = alpha_scale * rel_abund
            subj.reads[t] = dirichlet_multinomial(alpha, num_reads)

            # Make biomass
            triplicates = np.exp(
                np.log(total_mass[tidx]) +
                qpcr_noise_scale * pylab.random.normal.sample(size=3)
            )
            subj.qpcr[t] = pylab.qPCRdata(cfus=triplicates, mass=1., dilution_factor=1.)
    return study


def dirichlet_multinomial(alpha: np.ndarray, n: int) -> np.ndarray:
    return np.random.multinomial(n, np.random.dirichlet(alpha))


def main():
    args = parse_args()
    growth_rates, interactions, interaction_indicators, taxa_names, initial_cond_mean, initial_cond_std = parse_glv_params(Path(args.input_glv_params))
    time_points = parse_time_points(args.time_points_file)
    seed = args.seed

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    taxa = TaxaSet()
    for taxa_name in taxa_names:
        taxa.add_taxon(taxa_name)

    synthetic = make_synthetic('cdiff_mdsine_bvs', taxa, growth_rates, interactions, interaction_indicators, seed=seed)

    # Make subject names
    synthetic.set_subjects([f'subj_{i}' for i in range(args.num_subjects)] + ['holdout'])
    synthetic.set_timepoints(time_points)

    # Generate the trajectories.
    raw_trajs = synthetic.generate_trajectories(
        dt=args.sim_dt,
        init_dist=variables.Normal(initial_cond_mean, initial_cond_std),
        processvar=model.MultiplicativeGlobal(args.process_var)
    )

    # Plot the trajectories.
    for subj in synthetic.subjs:
        fig, ax = plt.subplots(figsize=(10, 8))
        trajs = raw_trajs[subj]['X']  # (n_taxa) x (n_times)
        times = raw_trajs[subj]['times']
        for taxa_traj in trajs:
            ax.plot(times, taxa_traj, marker=None)

        traj_plot_path = out_dir / f'{subj}.pdf'
        plt.savefig(traj_plot_path)
        plt.close(fig)
        print(f"Saved trajectories to {traj_plot_path}")

    # Simulate noise levels.
    noise_levels = {
        'low': args.low_noise,
        'medium': args.medium_noise,
        'high': args.high_noise
    }

    for noise_level_name, noise_level in noise_levels.items():
        print(f"Simulating QPcr noise level {noise_level_name}: {noise_level}, and reads noise model `DMD`.")

        # ======== Simulate noise using DMD to test robustness.
        study = simulate_reads_dmd(
            synth=synthetic,
            study_name=f'simulated-{noise_level_name}',
            alpha_scale=args.dmd_alpha_scale,
            num_reads=args.read_depth,
            qpcr_noise_scale=noise_level
        )
        study.perturbations = []

        # ======== Simulate using NegBin noise model.
        # study = synthetic.simulateMeasurementNoise(
        #     a0=args.negbin_a0,
        #     a1=args.negbin_a1,
        #     qpcr_noise_scale=noise_level,
        #     approx_read_depth=args.read_depth,
        #     name=f'simulated-{noise_level_name}'
        # )

        # Save to disk.
        holdout_study = study.pop_subject('holdout', name=f'holdout-{noise_level_name}')
        study.perturbations = None
        holdout_study.perturbations = None

        holdout_study.save(str(out_dir / f'holdout_{noise_level_name}.pkl'))
        study.save(str(out_dir / f'subjset_{noise_level_name}.pkl'))


if __name__ == "__main__":
    main()
