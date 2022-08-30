"""
Analyzes a keystonness-like array of eigenvalue decompositions.
"""

import argparse
from pathlib import Path
from typing import Iterator, Tuple

import numpy as np
import scipy.linalg

import mdsine2 as md2
from mdsine2 import TaxaSet
from mdsine2.names import STRNAMES
from mdsine2.base.cluster import _Cluster


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--fixed-cluster-mcmc-path', '-f', type=str, dest='mcmc_path',
                        required=True,
                        help='Path of saved MDSINE2.BaseMCMC chain (fixed-clustering inference)')
    parser.add_argument('--output', '-o', type=str, dest='out_dir',
                        required=True,
                        help='Directory to save desired output eigenvalues (.npz format).')
    return parser.parse_args()


def sub_interaction_matrix(interactions: np.ndarray, taxaset: TaxaSet, excluded_cluster: _Cluster) -> Iterator[np.ndarray]:
    excluded_indices = set(oidx for oidx in excluded_cluster.members)
    remaining_indices = [i for i in range(len(taxaset)) if i not in excluded_indices]
    rows, cols = np.ix_(remaining_indices, remaining_indices)

    for interaction_sample in interactions:
        yield interaction_sample[rows, cols]


def load_interactions(mcmc: md2.BaseMCMC):
    self_interactions = mcmc.graph[STRNAMES.SELF_INTERACTION_VALUE].get_trace_from_disk(section="posterior")
    interactions = mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk(section="posterior")
    interactions[np.isnan(interactions)] = 0
    self_interactions = -np.absolute(self_interactions)
    for i in range(self_interactions.shape[1]):
        interactions[:, i, i] = self_interactions[:, i]
    return interactions


def compute_eigenvalues(mcmc: md2.BaseMCMC, taxa: TaxaSet) -> Iterator[Tuple[int, int, np.ndarray]]:
    interactions = load_interactions(mcmc)
    for cluster_idx, cluster in enumerate(mcmc.graph[STRNAMES.CLUSTERING_OBJ]):
        eigs = []
        for interaction_matrix in sub_interaction_matrix(interactions, taxa, excluded_cluster=cluster):
            eig = scipy.linalg.eigvals(interaction_matrix)
            eigs.append(eig)
        yield cluster_idx, cluster.id, np.array(eigs)


def main():
    args = parse_args()

    mcmc_path = Path(args.mcmc_path)
    mcmc = md2.BaseMCMC.load(str(mcmc_path))
    study = md2.BaseMCMC.load(str(mcmc_path.parent / "subjset.pkl"))
    taxa = study.taxa

    out_dir = Path(args.out_dir)
    print(f"Output directory: {out_dir}")
    for cluster_idx, cluster_id, sample_eigs in compute_eigenvalues(mcmc, taxa):
        print(f"Computed eigenvalues by excluding cluster IDX:{cluster_idx} ID:{cluster_id}")
        out_path = out_dir / f'eigenvalues_exclude_cluster_{cluster_idx}-{cluster_id}.npz'
        np.save(str(out_path), sample_eigs)


if __name__ == "__main__":
    main()
