import argparse
import os
from pathlib import Path

import mdsine2 as md2
from mdsine2.names import STRNAMES


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input', '-i', type=str, dest='input',
        required=True,
        help='This is the dataset to do inference with.'
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
        '--log-every', type=int, default=100,
        required=False,
        help='<Optional> Tells the inference loop to print debug messages every k iterations.'
    )

    parser.add_argument(
        '--time_mask', type=str,
        required=False,
        help='<Optional> Specify spike-in timepoint for each taxa.'
    )
    return parser.parse_args()


def main():
    args = parse_args()
    basepath = Path(args.basepath)

    # 1) load dataset
    study = md2.Study.load(args.input)
    if study.perturbations is not None and len(study.perturbations) == 0:
        study.perturbations = None
    md2.seed(args.seed)

    # 2) Load the model parameters
    output_dir = basepath / study.name
    output_dir.mkdir(exist_ok=True, parents=True)

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

    # 3) Begin inference
    params = md2.config.MDSINE2ModelConfig(
        basepath=str(basepath), seed=args.seed,
        burnin=args.burnin, n_samples=args.n_samples, negbin_a1=a1,
        negbin_a0=a0, checkpoint=args.checkpoint)
    params.LEARN[STRNAMES.CLUSTERING] = False
    params.INITIALIZATION_KWARGS[STRNAMES.CLUSTERING]['value_option'] = 'no-clusters'

    if args.time_mask is not None and len(args.time_mask) > 0:
        params.INITIALIZATION_KWARGS[STRNAMES.ZERO_INFLATION]['value_option'] = "custom"
        params.ZERO_INFLATION_TRANSITION_POLICY = 'ignore'
        params.ZERO_INFLATION_DATA_PATH = Path(args.time_mask)

    # Run with multiprocessing if necessary
    if args.mp:
        params.MP_FILTERING = 'full'
        params.MP_CLUSTERING = 'full-4'

    # Set the sparsities
    params.INITIALIZATION_KWARGS[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB]['hyperparam_option'] = \
        args.interaction_prior
    params.INITIALIZATION_KWARGS[STRNAMES.PERT_INDICATOR_PROB]['hyperparam_option'] = \
        args.perturbation_prior

    mcmc = md2.initialize_graph(params=params, graph_name=study.name, subjset=study)
    mdata_fname = os.path.join(params.MODEL_PATH, 'metadata.txt')
    params.make_metadata_file(fname=mdata_fname)
    md2.run_graph(mcmc, crash_if_error=True, log_every=args.log_every)


if __name__ == "__main__":
    main()
