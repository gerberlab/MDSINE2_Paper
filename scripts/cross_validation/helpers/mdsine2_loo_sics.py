import os
import argparse

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
        '--nomodules', action='store_true', dest='nomodules',
        help='If flag is provided, then run inference without learning modules.'
    )
    parser.add_argument(
        '--exclude-subject', '-e', type=str, dest='exclude_subject',
        required=True,
        help='The name of the subject to exclude.'
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
        '--interaction-ind-prior-a', type=float, dest='interaction_prior_a',
        required=False, help='The `a` parameter of the interaction indicator prior.', default=None
    )
    parser.add_argument(
        '--interaction-ind-prior-b', type=float, dest='interaction_prior_b',
        required=False, help='The `b` parameter of the interaction indicator prior.', default=None
    )
    # =========== DEFAULT VALUES LISTED BELOW
    parser.add_argument(
        '--interaction-mean-loc', type=float, dest='interaction_mean_loc',
        required=False, help='The loc parameter for the interaction strength prior mean.', default=0.0
    )
    parser.add_argument(
        '--interaction-var-dof', type=float, dest='interaction_var_dof',
        required=False, help='The dof parameter for the interaction strength prior var.', default=2.01
    )
    parser.add_argument(
        '--interaction-var-rescale', type=float, dest='interaction_var_rescale',
        required=False,
        help='Controls the scale parameter for the interaction strength prior var, using the formula [SCALE]*E^2',
        default=1e-8
    )
    parser.add_argument(
        '--growth-var-rescale', type=float, dest='growth_var_rescale',
        required=False,
        help='Controls the initialization of the growth variance (Note: its not learned by default)'
             'Uses the fomrmula 1e4*[SCALE]*<default_value>',
        default=1.0
    )
    parser.add_argument(
        '--si-var-rescale', type=float, dest='si_var_rescale',
        required=False,
        help='Controls the initialization of the si variance (Note: its not learned by default)'
             'Uses the fomrmula 1e4*[SCALE]*<default_value>',
        default=1.0
    )
    parser.add_argument(
        '--pert-var-mean', type=float, dest='pert_var_mean',
        required=False,
        help='Controls the mean of the perturbation variance. DOF=2.01, scale=<mean>*(dof-2)/dof.',
        default=1e4
    )

    parser.add_argument(
        '--log-every', type=int, default=100,
        required=False,
        help='<Optional> Tells the inference loop to print debug messages every k iterations.'
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # 1) load dataset
    print('Loading dataset {}'.format(args.input))
    study = md2.Study.load(args.input)
    md2.seed(args.seed)

    # 1.5) Pop out the subject.
    print("Removing subject: ")
    from pathlib import Path
    base_dir = Path(args.basepath) / study.name
    base_dir.mkdir(exist_ok=True, parents=True)

    try:
        heldout_study = study.pop_subject(args.exclude_subject, f'cv-exclude-{args.exclude_subject}')
        heldout_study.save(str(base_dir / f'cv-exclude-{args.exclude_subject}.pkl'))
    except ValueError as e:
        print("Encountered ValueError. Available subjects: {}".format(
            list(study._subjects.keys())
        ))
        exit(1)

    # 2) Load the model parameters
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

    print('Setting a0 = {:.4E}, a1 = {:.4E}'.format(a0, a1))

    # 3) Begin inference
    params = md2.config.MDSINE2ModelConfig(
        basepath=str(base_dir), seed=args.seed,
        burnin=args.burnin, n_samples=args.n_samples, negbin_a1=a1,
        negbin_a0=a0, checkpoint=args.checkpoint)

    if args.nomodules:
        params.LEARN[STRNAMES.CLUSTERING] = False
        params.LEARN[STRNAMES.CONCENTRATION] = False
        params.INITIALIZATION_KWARGS[STRNAMES.CLUSTERING]['value_option'] = 'no-clusters'
        params.INITIALIZATION_KWARGS[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB]['N'] = 'fixed-clustering'
        params.INITIALIZATION_KWARGS[STRNAMES.PERT_INDICATOR_PROB]['N'] = 'fixed-clustering'

    # Set the sparsities
    if args.interaction_prior is None:
        raise ValueError('Must specify `--interaction-ind-prior`')
    else:
        params.INITIALIZATION_KWARGS[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB]['hyperparam_option'] = args.interaction_prior
    params.INITIALIZATION_KWARGS[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB]['a'] = args.interaction_prior_a
    params.INITIALIZATION_KWARGS[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB]['b'] = args.interaction_prior_b

    if args.perturbation_prior is None:
        raise ValueError('Must specify `--perturbation-ind-prior`')
    else:
        params.INITIALIZATION_KWARGS[STRNAMES.PERT_INDICATOR_PROB]['hyperparam_option'] = args.perturbation_prior

    # Set the interaction strength priors
    if args.interaction_mean_loc != 0.0:
        params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_MEAN_INTERACTIONS]['loc_option'] = 'manual'
        params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_MEAN_INTERACTIONS]['loc'] = args.interaction_mean_loc

    if args.interaction_var_dof != None:
        params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_VAR_INTERACTIONS]['dof_option'] = 'manual'
        params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_VAR_INTERACTIONS]['dof'] = args.interaction_var_dof

    params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_VAR_INTERACTIONS]['scale_option'] = 'inflated-median'
    params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_VAR_INTERACTIONS]['inflation_factor'] = args.interaction_var_rescale

    params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_VAR_SELF_INTERACTIONS]['inflation_factor'] = 1e4 * args.si_var_rescale
    params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_VAR_GROWTH]['inflation_factor'] = 1e4 * args.growth_var_rescale

    # perturbation variance
    params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_VAR_PERT]['dof_option'] = 'manual'
    params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_VAR_PERT]['dof'] = 2.01
    params.INITIALIZATION_KWARGS[STRNAMES.PRIOR_VAR_PERT]['target_mean'] = args.pert_var_mean


    # Change the cluster initialization to no clustering if there are less than 30 clusters
    if len(study.taxa) <= 30:
        print(
            'Since there is less than 30 taxa, we set the initialization of the clustering to `no-clusters`')
        params.INITIALIZATION_KWARGS[STRNAMES.CLUSTERING]['value_option'] = 'no-clusters'

    mcmc = md2.initialize_graph(params=params, graph_name=study.name, subjset=study)
    mdata_fname = os.path.join(params.MODEL_PATH, 'metadata.txt')
    params.make_metadata_file(fname=mdata_fname)

    _ = md2.run_graph(mcmc, crash_if_error=True, log_every=args.log_every)


if __name__ == "__main__":
    main()
