"""
Analyzes a keystonness-like array of eigenvalue decompositions.
"""

import argparse
from pathlib import Path
from typing import Iterator, Tuple, Set, List

import numpy as np
import scipy.linalg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs-dir', '-i', type=str, dest='inputs_dir', required=True)
    parser.add_argument('--out-dir', '-o', type=str, dest='out_dir', required=True)
    return parser.parse_args()


def sub_interaction_matrix(interactions: np.ndarray, excluded_module: Set[int]) -> Iterator[np.ndarray]:
    n_taxa = interactions.shape[1]
    remaining_indices = [i for i in range(n_taxa) if i not in excluded_module]
    rows, cols = np.ix_(remaining_indices, remaining_indices)

    for interaction_sample in interactions:
        yield interaction_sample[rows, cols]


def compute_eigenvalues(interactions: np.ndarray, modules: List[Set[int]]) -> Iterator[Tuple[int, int, np.ndarray]]:
    for module_idx, module in enumerate(modules):
        eigs = []
        for interaction_matrix in sub_interaction_matrix(interactions, excluded_module=module):
            eig = scipy.linalg.eigvals(interaction_matrix)
            eigs.append(eig)
        yield module_idx, np.array(eigs)


def load_consensus_clustering(inputs_dir: Path) -> List[Set[int]]:
    agglom = np.load(str(inputs_dir / "agglomeration.npy"))

    clusters = []
    for cidx in range(np.max(agglom) + 1):  # Make sure to do the (+1) to count the last module.
        cluster = set(int(x) for x in np.where(agglom == cidx)[0])
        clusters.append(cluster)
    return clusters


def main():
    args = parse_args()
    inputs_dir = Path(args.inputs_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    modules: List[Set[int]] = load_consensus_clustering(inputs_dir)
    interactions = np.load(str(inputs_dir / "interactions.npy"))
    print(f"Output directory: {out_dir}")
    for cluster_idx, sample_eigs in compute_eigenvalues(interactions, modules):
        print(f"Computed eigenvalues by excluding cluster IDX:{cluster_idx}")
        out_path = out_dir / f'eigenvalues_exclude_cluster_{cluster_idx}.npy'
        np.save(str(out_path), sample_eigs)


if __name__ == "__main__":
    main()
