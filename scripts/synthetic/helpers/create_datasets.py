"""
Python script for generating semisynthetic samples for a given seed + noise level.
Takes as input MDSINE1's BVS sample matrix file
"""
from dataclasses import dataclass
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
    parser.add_argument('-p', '--perturbations_file', type=str, required=False, default=".",
                        help='<Optional> Specify a .json perturbations, representing a list of perturbation objects'
                             '{`name`: <str>, `strength`: <array[float]>, `start`: float, `end`: float}')
    parser.add_argument('-n', '--num_subjects', type=int, required=True,
                        help='<Required> The number of subjecs to simulate to lump into a single cohort.')
    parser.add_argument('-o', '--out_dir', type=str, required=True,
                        help='<Required> The directory to output the sampled subject to.')
    parser.add_argument('-qt', '--num_qpcr_triplicates', type=int, required=True,
                        help='<Required> The number of qPCR triplicates to sample.')
    parser.add_argument('-s', '--seed', type=int, required=True,
                        help='<Required> The seed to use for random sampling.')

    # Optional parameters
    parser.add_argument('-pv', '--process_var', type=float, required=False, default=0.01)
    parser.add_argument('-dt', '--sim_dt', type=float, required=False, default=0.01)
    parser.add_argument('-a', '--dmd_alpha_scale', type=float, required=False, default=286,
                        help='DMD dispersion parameter estimated from mean estimates at all time-points from '
                             'C. diff data (Default: 286, carry-over from MDSINE1)')
    parser.add_argument('-q', '--qpcr_noise', type=float, required=False, default=0.01)
    parser.add_argument('-r', '--read_depth', type=int, required=False, default=50000)
    parser.add_argument('-in', '--intervene_day', dest='intervene_day', type=float, required=False, default=0.0)
    parser.add_argument('-min', '--initial_min_value', type=float, required=False, default=100.0)

    # VARIANCE SCALING (noise levels)
    parser.add_argument('--low_noise', type=float, required=False, default=0.5)
    parser.add_argument('--medium_noise', type=float, required=False, default=1.0)
    parser.add_argument('--high_noise', type=float, required=False, default=2.0)
    return parser.parse_args()


@dataclass
class ParsedPerturbation(object):
    name: str
    strength: np.ndarray
    start: float
    end: float


def make_synthetic(
        taxa: TaxaSet,
        growth_rate_values: np.ndarray,
        interaction_values: np.ndarray,
        interaction_indicators: np.ndarray,
        perturbations: List[ParsedPerturbation],
        seed: int
) -> Synthetic:
    syn = Synthetic(name='SyntheticModel', seed=seed)
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
    syn.model.perturbations = [pert.strength for pert in perturbations]
    syn.model.perturbation_ends = [pert.end for pert in perturbations]
    syn.model.perturbation_starts = [pert.start for pert in perturbations]
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


def simulate_reads_dmd(synth: Synthetic, study_name: str, perts: List[ParsedPerturbation], alpha_scale: float, num_reads: int, qpcr_noise_scale: float) -> Study:
    # Make the study object
    study = Study(taxa=synth.taxa, name=study_name)
    for subjname in synth.subjs:
        study.add_subject(name=subjname)

    # Add times for each subject
    for subj in study:
        subj.times = synth.times

    # Add perts for each subject
    if synth.perturbations is not None and len(synth.perturbations) > 0:
        study.perturbations = Perturbations()
        for pert in perts:
            study.perturbations.append(BasePerturbation(
                name=pert.name,
                starts={subj.name: pert.start for subj in study},
                ends={subj.name: pert.end for subj in study}
            ))

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
            subj.qpcr[t] = qPCRdata(cfus=triplicates, mass=1., dilution_factor=1.)
    return study


def simulate_replicates(synth: Synthetic,
                        subj_idx: int,
                        subj_timepoints: List[int],  # indices
                        num_replicates: int,
                        taxa: TaxaSet,
                        study_name: str,
                        alpha_scale: float,
                        num_reads: int,
                        qpcr_noise_scale: float) -> Study:
    # Make the study object
    replicate_study = Study(taxa=taxa, name=study_name)
    for t in subj_timepoints:
        replicate_study.add_subject(name=f'M{subj_idx}-T{t}')

    # Add times for each subject
    replicate_times = np.arange(0., num_replicates, 1.0)
    for replicate_subj in replicate_study:
        replicate_subj.times = replicate_times

    for replicate_subj, subj_t in zip(replicate_study, subj_timepoints):
        # Extract timepoint abundances
        conc = synth._data[synth.subjs[subj_idx]][:, subj_t]
        total_mass = np.sum(conc)
        rel_mass = conc / total_mass
        for repl_t in replicate_times:
            # Make the reads
            alpha = alpha_scale * rel_mass
            replicate_subj.reads[repl_t] = dirichlet_multinomial(alpha, num_reads)

            # Make biomass
            triplicates = np.exp(
                np.log(total_mass) +
                qpcr_noise_scale * pylab.random.normal.sample(size=3)
            )
            replicate_subj.qpcr[repl_t] = qPCRdata(cfus=triplicates, mass=1., dilution_factor=1.)
    return replicate_study


def dirichlet_multinomial(alpha: np.ndarray, n: int) -> np.ndarray:
    support = alpha[alpha > 0]
    x = np.zeros(len(alpha), dtype=int)
    x[alpha > 0] = np.random.multinomial(n, np.random.dirichlet(support))
    return x


def simulate_trajectories(synth: Synthetic,
                          init_dist: Variable,
                          taxa: TaxaSet,
                          initial_min_value: float,
                          dt: float,
                          processvar: model.MultiplicativeGlobal,
                          intervene_day: float = 0.0):
    raw_trajs = {}

    for subj in synth.subjs:
        print('Forward simulating {}'.format(subj))
        init_abund = init_dist.sample(size=len(taxa))
        init_abund[init_abund < initial_min_value] = initial_min_value

        if intervene_day == 0:
            total_n_days = synth.times[-1]
            d = pylab.integrate(
                dynamics=synth.model,
                initial_conditions=init_abund.reshape(-1, 1),
                dt=dt,
                n_days=total_n_days + dt,
                processvar=processvar,
                subsample=False
            )
        else:
            print(f"Simulating trajectory with intervention at day {intervene_day}")
            pathogen_abund = init_abund[0]
            init_abund[0] = 0.0

            total_n_days = synth.times[-1]
            d_pre = pylab.integrate(dynamics=synth.model, initial_conditions=init_abund.reshape(-1, 1),
                                    dt=dt, n_days=intervene_day, processvar=processvar,
                                    subsample=False)

            new_abund = d_pre['X'][:, -1]
            new_abund[0] = pathogen_abund
            d_post = pylab.integrate(dynamics=synth.model, initial_conditions=new_abund.reshape(-1, 1),
                                     dt=dt, n_days=total_n_days - intervene_day + dt, processvar=processvar,
                                     subsample=False)
            # Merge the two results
            d = {
                'X': np.concatenate([d_pre['X'], d_post['X']], axis=1),
                'times': np.concatenate([d_pre['times'], d_post['times'] + intervene_day], axis=0)
            }

        # Save the trajectories
        steps_per_day = int(np.ceil(total_n_days / dt) / total_n_days)
        idxs = [int(steps_per_day * t) for t in synth.times]
        X = d['X']
        synth._data[subj] = X[:, idxs]
        raw_trajs[subj] = d
    return raw_trajs


def parse_perturbations(pert_path: Path) -> List[ParsedPerturbation]:
    if not pert_path.exists() or pert_path.is_dir():
        return []

    import json
    with open(pert_path, "r") as pf:
        raw_objs = json.load(pf)
        return [
            ParsedPerturbation(
                name=obj['name'],
                strength=np.array(obj['strength']),
                start=obj['start'],
                end=obj['end']
            )
            for obj in raw_objs
        ]


def main():
    args = parse_args()
    growth_rates, interactions, interaction_indicators, taxa_names, initial_cond_mean, initial_cond_std = parse_glv_params(Path(args.input_glv_params))
    perturbations = parse_perturbations(Path(args.perturbations_file))
    time_points = parse_time_points(args.time_points_file)
    seed = args.seed

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    taxa = TaxaSet()
    for taxa_name in taxa_names:
        taxa.add_taxon(taxa_name)

    synthetic = make_synthetic(taxa, growth_rates, interactions, interaction_indicators, perturbations, seed=seed)

    # Make subject names
    synthetic.set_subjects([f'subj_{i}' for i in range(args.num_subjects)])
    synthetic.set_timepoints(time_points)

    if args.process_var > 0:
        process_var = model.MultiplicativeGlobal(args.process_var)
    else:
        process_var = None

    # Generate the trajectories.
    raw_trajs = simulate_trajectories(
        synth=synthetic,
        taxa=taxa,
        dt=args.sim_dt,
        init_dist=variables.Normal(initial_cond_mean, np.power(initial_cond_std, 2)),
        processvar=process_var,
        initial_min_value=args.initial_min_value,
        intervene_day=args.intervene_day
    )

    # Plot the trajectories.
    for subj in synthetic.subjs:
        fig, ax = plt.subplots(figsize=(10, 8))
        trajs = raw_trajs[subj]['X']  # (n_taxa) x (n_times)
        times = raw_trajs[subj]['times']
        for taxa_traj in trajs:
            ax.plot(times, taxa_traj, marker=None)
        ax.set_yscale('log')

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

    # Sample the data.
    for noise_level_name, variance_scaling in noise_levels.items():
        dmd_scale = args.dmd_alpha_scale / variance_scaling
        print(f"Simulating DMD noise level {noise_level_name}: {dmd_scale}.")

        # ======== Simulate noise using DMD to test robustness.
        study = simulate_reads_dmd(
            synth=synthetic,
            study_name=f'simulated-{noise_level_name}',
            perts=perturbations,
            alpha_scale=dmd_scale,
            num_reads=args.read_depth,
            qpcr_noise_scale=args.qpcr_noise
        )

        replicate_study = simulate_replicates(
            synth=synthetic,
            subj_idx=0,
            subj_timepoints=[13, 19, 25],
            num_replicates=args.num_qpcr_triplicates,
            taxa=taxa,
            study_name=f'replicate-{noise_level_name}',
            alpha_scale=dmd_scale,
            num_reads=args.read_depth,
            qpcr_noise_scale=args.qpcr_noise
        )

        pkl_dir = out_dir / f'reads_{args.read_depth}' / f'noise_{noise_level_name}'
        pkl_dir.mkdir(exist_ok=True, parents=True)

        study.save(str(pkl_dir / f'subjset.pkl'))
        print("Generated main subjset.")

        replicate_study.save(str(pkl_dir / f'replicate.pkl'))
        print("Generated qPCR replicates for NegBin fitting.")


if __name__ == "__main__":
    main()
