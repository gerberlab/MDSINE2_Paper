"""
Analyzes a keystonness-like array of eigenvalue decompositions.
"""

import argparse

import mdsine2 as md2
from mdsine2.names import STRNAMES
from mdsine2.pylab.cluster import _Cluster


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--fixed-cluster-mcmc-path', '-f', type=str, dest='mcmc_path',
                        required=True,
                        help='Path of saved MDSINE2.BaseMCMC chain (fixed-clustering inference)')
    parser.add_argument('--output', '-o', type=str, dest='out_path',
                        required=True,
                        help='Path to desired output eigenvalues (.npz format).')
    return parser.parse_args()


def sub_interaction_matrix(mcmc: md2.BaseMCMC, excluded_cluster: _Cluster) -> np.ndarray:
    pass


def main():
    args = parse_args()
    mcmc = md2.BaseMCMC.load(args.mcmc_path)

    for cluster_idx, cluster in enumerate(mcmc.graph[STRNAMES.CLUSTERING_OBJ]):
        interaction_matrix = sub_interaction_matrix(mcmc, excluded_cluster=cluster)

        # Compute eigenvalues



if __name__ == "__main__":
    main()
