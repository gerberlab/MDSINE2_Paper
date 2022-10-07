import os
import time
import argparse
from pathlib import Path
from typing import List

import numpy as np
import scipy.stats
from sklearn.cluster import AgglomerativeClustering

import mdsine2 as md2
from mdsine2.logger import logger
from mdsine2.names import STRNAMES


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input', '-i', type=str, dest='input',
        required=True,
        help='This is the dataset to do inference with.'
    )
    parser.add_argument(
        '--multiseed-dir', '-m', type=str, dest='multiseed_dir',
        required=True,
        help='The directory containing the .npy files for the combined seeds. Requires the '
    )
    parser.add_argument(
        '--negbin', type=str, dest='negbin', nargs='+',
        required=True,
        help='If there is a single argument, then this is the MCMC object that was run to ' \
             'learn a0 and a1. If there are two arguments passed, these are the a0 and a1 ' \
             'of the negative binomial dispersion parameters. Example: ' \
             '--negbin /path/to/negbin/mcmc.pkl. Example: ' \
             '--negbin 0.0025 0.025'
    )
    parser.add_argument(
        '--seed', '-s', type=int, dest='seed',
        required=True,
        help='This is the seed to initialize the inference with'
    )
    parser.add_argument(
        '--burnin', '-nb', type=int, dest='burnin',
        required=True,
        help='How many burn-in Gibb steps for Markov Chain Monte Carlo (MCMC)'
    )
    parser.add_argument(
        '--n-samples', '-ns', type=int, dest='n_samples',
        required=True,
        help='Total number Gibb steps to perform during MCMC inference'
    )
    parser.add_argument(
        '--checkpoint', '-c', type=int, dest='checkpoint',
        required=True,
        help='How often to write the posterior to disk. Note that `--burnin` and ' \
             '`--n-samples` must be a multiple of `--checkpoint` (e.g. checkpoint = 100, ' \
             'n_samples = 600, burnin = 300)'
    )
    parser.add_argument(
        '--basepath', '--output-basepath', '-b', type=str, dest='basepath',
        required=True,
        help='This is folder to save the output of inference'
    )
    parser.add_argument(
        '--multiprocessing', '-mp', type=int, dest='mp',
        help='If 1, run the inference with multiprocessing. Else run on a single process',
        default=0
    )
    parser.add_argument(
        '--rename-study', type=str, dest='rename_study',
        required=False, default=None,
        help='Specify the name of the study to set'
    )
    parser.add_argument(
        '--interaction-ind-prior', '-ip', type=str, dest='interaction_prior',
        required=True,
        help='Prior of the indicator of the interactions.'
    )
    parser.add_argument(
        '--perturbation-ind-prior', '-pp', type=str, dest='perturbation_prior',
        required=True,
        help='Prior of the indicator of the perturbations'
    )

    parser.add_argument(
        '--log-every', type=int, default=100, dest='log_every',
        required=False,
        help='<Optional> Tells the inference loop to print debug messages every k iterations.'
    )

    parser.add_argument(
        '--benchmark', action='store_true', dest='benchmark',
        help='If flag is set, then logs (at INFO level) the update() runtime of each component at the end.'
    )
    return parser.parse_args()


def agglomerate_from_cocluster(multiseed_dir: Path) -> List[List[int]]:
    A = 1 - np.load(str(multiseed_dir / "coclusters.npy"))
    n = scipy.stats.mode(np.load(str(multiseed_dir / "n_clusters.npy")))[0][0]

    linkage = 'average'
    c = AgglomerativeClustering(
        n_clusters=n,
        affinity='precomputed',
        linkage=linkage
    )

    agglom = c.fit_predict(A)
    clusters = []
    for cidx in range(np.max(agglom) + 1):  # Make sure to do the (+1) to count the last module.
        cluster = list(np.where(agglom == cidx)[0])
        clusters.append([int(oidx) for oidx in cluster])
    return clusters


def main():
    args = parse_args()

    # 1) load dataset
    logger.info('Loading dataset {}'.format(args.input))
    study = md2.Study.load(args.input)
    if args.rename_study is not None:
        if args.rename_study.lower() != 'none':
            study.name = args.rename_study
    md2.seed(args.seed)

    # 2) Load the model parameters
    os.makedirs(args.basepath, exist_ok=True)
    basepath = os.path.join(args.basepath, study.name)
    os.makedirs(basepath, exist_ok=True)

    # Load the negative binomial parameters
    if len(args.negbin) == 1:
        negbin = md2.BaseMCMC.load(args.negbin[0])
        a0 = md2.summary(negbin.graph[STRNAMES.NEGBIN_A0])['mean']
        a1 = md2.summary(negbin.graph[STRNAMES.NEGBIN_A1])['mean']

    elif len(args.negbin) == 2:
        a0 = float(args.negbin[0])
        a1 = float(args.negbin[1])
    else:
        raise ValueError('Argument `negbin`: there must be only one or two arguments.')

    logger.info('Setting a0 = {:.4E}, a1 = {:.4E}'.format(a0, a1))

    # 2.5) Load agglomerated modules
    agglomerated_modules = agglomerate_from_cocluster(Path(args.multiseed_dir))

    # 3) Begin inference
    params = md2.config.MDSINE2ModelConfig(
        basepath=basepath, seed=args.seed,
        burnin=args.burnin, n_samples=args.n_samples, negbin_a1=a1,
        negbin_a0=a0, checkpoint=args.checkpoint)
    # Run with multiprocessing if necessary
    if args.mp:
        params.MP_FILTERING = 'full'
        params.MP_CLUSTERING = 'full-4'

    # Load fixed clustering structure
    params.LEARN[STRNAMES.CLUSTERING] = False
    params.LEARN[STRNAMES.CONCENTRATION] = False
    params.INITIALIZATION_KWARGS[STRNAMES.CLUSTERING]['value_option'] = 'manual'
    params.INITIALIZATION_KWARGS[STRNAMES.CLUSTERING]['value'] = agglomerated_modules
    params.INITIALIZATION_KWARGS[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB]['N'] = 'fixed-clustering'
    params.INITIALIZATION_KWARGS[STRNAMES.PERT_INDICATOR_PROB]['N'] = 'fixed-clustering'

    # Set the sparsities
    params.INITIALIZATION_KWARGS[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB]['hyperparam_option'] = \
        args.interaction_prior
    params.INITIALIZATION_KWARGS[STRNAMES.PERT_INDICATOR_PROB]['hyperparam_option'] = \
        args.perturbation_prior

    # Change the cluster initialization to no clustering if there are less than 30 clusters
    if len(study.taxa) <= 30:
        logger.info(
            'Since there is less than 30 taxa, we set the initialization of the clustering to `no-clusters`')
        params.INITIALIZATION_KWARGS[STRNAMES.CLUSTERING]['value_option'] = 'no-clusters'

    mcmc = md2.initialize_graph(params=params, graph_name=study.name, subjset=study)
    mdata_fname = os.path.join(params.MODEL_PATH, 'metadata.txt')
    params.make_metadata_file(fname=mdata_fname)

    start_time = time.time()
    mcmc = md2.run_graph(mcmc, crash_if_error=True, log_every=args.log_every, benchmarking=args.benchmark)

    # Record how much time inference took
    t = time.time() - start_time
    t = t / 3600  # Convert to hours

    f = open(mdata_fname, 'a')
    f.write('\n\nTime for inference: {} hours'.format(t))
    f.close()


if __name__ == "__main__":
    main()
