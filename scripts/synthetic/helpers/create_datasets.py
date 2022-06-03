"""
Python script for generating semisynthetic samples for a given seed + noise level.
Takes as input MDSINE1's BVS sample matrix file
"""
from typing import Tuple, List
from pathlib import Path
import argparse

import numpy as np
from mdsine2 import *
from mdsine2.names import STRNAMES


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_glv_params', type=str, required=True,
                        help='<Required> The path to a .npz file containing `growth_rates` and `interactions` arrays')
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
    parser.add_argument('-a0', '--negbin_a0', type=float, required=False, default=1e-10)
    parser.add_argument('-a1', '--negbin_a1', type=float, required=False, default=0.05)
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
    initial_cond_var = params['initial_std'] * np.ones(len(growth), dtype=float)
    indicators = (interactions != 0.0)
    taxa_names = [f'TAXA_{i+1}' for i in range(len(growth))]
    return growth, interactions, indicators, taxa_names, initial_cond_mean, initial_cond_var


def parse_time_points(time_points_path: Path) -> np.ndarray:
    time_points = []
    with open(time_points_path, 'rt') as f:
        for line in f:
            t = float(line.rstrip())
            time_points.append(t)
    return np.array(time_points, dtype=float)


def main():
    args = parse_args()
    growth_rates, interactions, interaction_indicators, taxa_names, initial_cond_mean, initial_cond_var = parse_glv_params(Path(args.input_glv_params))
    time_points = parse_time_points(args.time_points_file)
    seed = args.seed

    taxa = TaxaSet()
    for taxa_name in taxa_names:
        taxa.add_taxon(taxa_name)

    synthetic = make_synthetic('cdiff_mdsine_bvs', taxa, growth_rates, interactions, interaction_indicators, seed=seed)

    # Make subject names
    synthetic.set_subjects([f'subj_{i}' for i in range(args.num_subjects)])
    synthetic.set_timepoints(time_points)

    # Generate the trajectories.
    synthetic.generate_trajectories(
        dt=args.sim_dt,
        init_dist=variables.Normal(initial_cond_mean, initial_cond_var),
        processvar=model.MultiplicativeGlobal(args.process_var)
    )

    noise_levels = {
        'low': args.low_noise,
        'medium': args.medium_noise,
        'high': args.high_noise
    }

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)
    for noise_level_name, noise_level in noise_levels.items():
        # Simulate noise.
        study = synthetic.simulateMeasurementNoise(
            a0=args.negbin_a0,
            a1=args.negbin_a1,
            qpcr_noise_scale=noise_level,
            approx_read_depth=args.read_depth,
            name=f'cdiff-semisynth-noise-{noise_level_name}'
        )
        study.save(str(out_dir / f'subjset_{noise_level_name}.pkl'))


if __name__ == "__main__":
    main()
