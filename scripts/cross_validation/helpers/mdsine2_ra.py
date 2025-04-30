"""
Implementation of MDSINE2 without time-series qPCR measurements.
"""

import os
from typing import Dict, Union, Tuple, List

import numpy as np
from mdsine2 import (
    config,
    pylab as pl,
    design_matrices,
    Study, BaseMCMC, Subject,
    Clustering, ClusterPerturbationEffect, posterior,
    qPCRdata
)
from mdsine2.names import STRNAMES
from mdsine2.logger import logger
from mdsine2.posterior import negbin_loglikelihood_MH_condensed_not_fast, negbin_loglikelihood_MH_condensed


class FilteringWithoutQPCR(posterior.FilteringLogMP):
    def __init__(self, biomass_mean: float, biomass_cv: float, *args, **kwargs):
        """
        "cv": "coefficient of variation" = stdev / mean
        """
        self.biomass_mean = biomass_mean
        self.biomass_cv = biomass_cv
        super().__init__(*args, **kwargs)

    def initialize(self, x_value_option: str, a0: Union[float, str], a1: Union[float, str],
        v1: Union[float, int], v2: Union[float, int], essential_timepoints: Union[np.ndarray, str],
        tune: Tuple[int, int], proposal_init_scale: Union[float, int],
        intermediate_step: Union[Tuple[str, Tuple], np.ndarray, List, type(None)],
        intermediate_interpolation: str=None, delay: int=0, bandwidth: Union[float, int]=None,
        window: int=None, target_acceptance_rate: Union[str, float]=0.44,
        calculate_qpcr_loglik: bool=True):
        '''
        Copy-pasted from posterior.py implementation
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay
        self._there_are_perturbations = self.G.perturbations is not None

        # Set the hyperparameters
        if not pl.isfloat(target_acceptance_rate):
            raise TypeError('`target_acceptance_rate` must be a float'.format(
                type(target_acceptance_rate)))
        if target_acceptance_rate < 0 or target_acceptance_rate > 1:
            raise ValueError('`target_acceptance_rate` ({}) must be in (0,1)'.format(
                target_acceptance_rate))
        if not pl.istuple(tune):
            raise TypeError('`tune` ({}) must be a tuple'.format(type(tune)))
        if len(tune) != 2:
            raise ValueError('`tune` ({}) must have 2 elements'.format(len(tune)))
        if not pl.isint(tune[0]):
            raise TypeError('`tune` ({}) 1st parameter must be an int'.format(type(tune[0])))
        if tune[0] < 0:
            raise ValueError('`tune` ({}) 1st parameter must be > 0'.format(tune[0]))
        if not pl.isint(tune[1]):
            raise TypeError('`tune` ({}) 2nd parameter must be an int'.format(type(tune[1])))
        if tune[1] < 0:
            raise ValueError('`tune` ({}) 2nd parameter must be > 0'.format(tune[1]))

        if not pl.isnumeric(a0):
            raise TypeError('`a0` ({}) must be a numeric type'.format(type(a0)))
        elif a0 <= 0:
            raise ValueError('`a0` ({}) must be > 0'.format(a0))
        if not pl.isnumeric(a1):
            raise TypeError('`a1` ({}) must be a numeric type'.format(type(a1)))
        elif a1 <= 0:
            raise ValueError('`a1` ({}) must be > 0'.format(a1))

        if not pl.isnumeric(proposal_init_scale):
            raise TypeError('`proposal_init_scale` ({}) must be a numeric type (int, float)'.format(
                type(proposal_init_scale)))
        if proposal_init_scale < 0:
            raise ValueError('`proposal_init_scale` ({}) must be positive'.format(
                proposal_init_scale))

        self.tune = tune
        self.a0 = a0
        self.a1 = a1
        self.target_acceptance_rate = target_acceptance_rate
        self.proposal_init_scale = proposal_init_scale
        self.v1 = v1
        self.v2 = v2

        # Set the essential timepoints (check to see if there is any missing data)
        if essential_timepoints is not None:
            logger.info('Setting up the essential timepoints')
            if pl.isstr(essential_timepoints):
                if essential_timepoints in ['auto', 'union']:
                    essential_timepoints = set()
                    for ts in self.G.data.times:
                        essential_timepoints = essential_timepoints.union(set(list(ts)))
                    essential_timepoints = np.sort(list(essential_timepoints))
                else:
                    raise ValueError('`essential_timepoints` ({}) not recognized'.format(
                        essential_timepoints))
            elif not pl.isarray(essential_timepoints):
                raise TypeError('`essential_timepoints` ({}) must be a str or an array'.format(
                    type(essential_timepoints)))
            logger.info('Essential timepoints: {}'.format(essential_timepoints))
            self.G.data.set_timepoints(times=essential_timepoints, eps=None, reset_timepoints=True)
            self.x.reset_value_size()

        # Set the intermediate timepoints if necessary
        if intermediate_step is not None:
            # Set the intermediate timepoints in the data
            if not pl.istuple(intermediate_step):
                raise TypeError('`intermediate_step` ({}) must be a tuple'.format(
                    type(intermediate_step)))
            if len(intermediate_step) != 2:
                raise ValueError('`intermediate_step` ({}) must be length 2'.format(
                    len(intermediate_step)))
            f, args = intermediate_step
            if not pl.isstr(f):
                raise TypeError('intermediate_step type ({}) must be a str'.format(type(f)))
            if f == 'step':
                if not pl.istuple(args):
                    raise TypeError('`args` ({}) must be a tuple'.format(type(args)))
                if len(args) != 2:
                    raise TypeError('`args` ({}) must have 2 arguments'.format(len(args)))
                step, eps = args
                self.G.data.set_timepoints(timestep=step, eps=eps, reset_timepoints=False)
            elif f == 'preserve-density':
                if not pl.istuple(args):
                    raise TypeError('`args` ({}) must be a tuple'.format(type(args)))
                if len(args) != 2:
                    raise TypeError('`args` ({}) must have 2 arguments'.format(len(args)))
                n, eps = args
                if not pl.isint(n):
                    raise TypeError('`n` ({}) must be an int'.format(type(n)))

                # For each timepoint, add `n` intermediate timepoints
                for ridx in range(self.G.data.n_replicates):
                    times = []
                    for i in range(len(self.G.data.times[ridx])-1):
                        t0 = self.G.data.times[ridx][i]
                        t1 = self.G.data.times[ridx][i+1]
                        step = (t1-t0)/(n+1)
                        times = np.append(times, np.arange(t0,t1,step=step))
                    times = np.sort(np.unique(times))
                    # print('\n\ntimes to put in', times)
                    self.G.data.set_timepoints(times=times, eps=eps, ridx=ridx, reset_timepoints=False)
                    # print('times for ridx {}'.format(self.G.data.times[ridx]))
                    # print('len times', len(self.G.data.times[ridx]))
                    # print('data shape', self.G.data.data[ridx].shape)

                # sys.exit()
            elif f == 'manual':
                raise NotImplementedError('Not Implemented')
            else:
                raise ValueError('`intermediate_step type ({}) not recognized'.format(f))
            self.x.reset_value_size()

        if intermediate_interpolation is not None:
            if intermediate_interpolation in ['linear-interpolation', 'auto']:
                for ridx in range(self.G.data.n_replicates):
                    for tidx in range(self.G.data.n_timepoints_for_replicate[ridx]):
                        if tidx not in self.G.data.given_timeindices[ridx]:
                            # We need to interpolate this time point
                            # get the previous given and next given timepoint
                            prev_tidx = None
                            for ii in range(tidx-1,-1,-1):
                                if ii in self.G.data.given_timeindices[ridx]:
                                    prev_tidx = ii
                                    break
                            if prev_tidx is None:
                                # Set to the same as the closest forward timepoint then continue
                                next_idx = None
                                for ii in range(tidx+1, self.G.data.n_timepoints_for_replicate[ridx]):
                                    if ii in self.G.data.given_timeindices[ridx]:
                                        next_idx = ii
                                        break
                                self.G.data.data[ridx][:,tidx] = self.G.data.data[ridx][:,next_idx]
                                continue

                            next_tidx = None
                            for ii in range(tidx+1, self.G.data.n_timepoints_for_replicate[ridx]):
                                if ii in self.G.data.given_timeindices[ridx]:
                                    next_tidx = ii
                                    break
                            if next_tidx is None:
                                # Set to the previous timepoint then continue
                                self.G.data.data[ridx][:,tidx] = self.G.data.data[ridx][:,prev_tidx]
                                continue

                            # Interpolate from prev_tidx to next_tidx
                            x = self.G.data.times[ridx][tidx]
                            x0 = self.G.data.times[ridx][prev_tidx]
                            y0 = self.G.data.data[ridx][:,prev_tidx]
                            x1 = self.G.data.times[ridx][next_tidx]
                            y1 = self.G.data.data[ridx][:,next_tidx]
                            self.G.data.data[ridx][:,tidx] = y0 * (1-((x-x0)/(x1-x0))) + y1 * (1-((x1-x)/(x1-x0)))
            else:
                raise ValueError('`intermediate_interpolation` ({}) not recognized'.format(intermediate_interpolation))

        # Initialize the latent trajectory
        if not pl.isstr(x_value_option):
            raise TypeError('`x_value_option` ({}) is not a str'.format(type(x_value_option)))
        if x_value_option == 'coupling':
            self._init_coupling()
        elif x_value_option == 'moving-avg':
            if not pl.isnumeric(bandwidth):
                raise TypeError('`bandwidth` ({}) must be a numeric'.format(type(bandwidth)))
            if bandwidth <= 0:
                raise ValueError('`bandwidth` ({}) must be positive'.format(bandwidth))
            self.bandwidth = bandwidth
            self._init_moving_avg()
        elif x_value_option in ['loess', 'auto']:
            if window is None:
                raise TypeError('If `value_option` is loess, then `window` must be specified')
            if not pl.isint(window):
                raise TypeError('`window` ({}) must be an int'.format(type(window)))
            if window <= 0:
                raise ValueError('`window` ({}) must be > 0'.format(window))
            self.window = window
            self._init_loess()
        else:
            raise ValueError('`x_value_option` ({}) not recognized'.format(x_value_option))

        # Get necessary data and set the parallel objects
        if self._there_are_perturbations:
            pert_starts = []
            pert_ends = []
            for perturbation in self.G.perturbations:
                pert_starts.append(perturbation.starts)
                pert_ends.append(perturbation.ends)
        else:
            pert_starts = None
            pert_ends = None

        if self.mp is None:
            self.mp = 'debug'
        if not pl.isstr(self.mp):
            raise TypeError('`mp` ({}) must either be a string or None'.format(type(self.mp)))
        if self.mp == 'debug':
            self.pool = []
        elif self.mp == 'full':
            self.pool = pl.multiprocessing.PersistentPool(G=self.G, ptype='sadw')
            self.worker_pids = []
        else:
            raise ValueError('`mp` ({}) not recognized'.format(self.mp))

        for ridx, subj in enumerate(self.G.data.subjects):
            # Set up qPCR measurements and reads to send
            qpcr_log_measurements = {}
            for t in self.G.data.given_timepoints[ridx]:
                qpcr_log_measurements[t] = self.G.data.qpcr[ridx][t].log_data
            reads = self.G.data.subjects.iloc(ridx).reads

            worker = SubjectLogTrajectorySetWithoutQPCR(biomass_mean=self.biomass_mean, biomass_cv=self.biomass_cv)
            worker.initialize(
                zero_inflation_transition_policy=self.zero_inflation_transition_policy,
                times=self.G.data.times[ridx],
                qpcr_log_measurements=qpcr_log_measurements,
                reads=reads,
                there_are_intermediate_timepoints=True,
                there_are_perturbations=self._there_are_perturbations,
                pv_global=self.G[STRNAMES.PROCESSVAR].global_variance,
                x_prior_mean=np.log(1e7),
                x_prior_std=1e10,
                tune=tune[1],
                delay=delay,
                end_iter=tune[0],
                proposal_init_scale=proposal_init_scale,
                a0=a0,
                a1=a1,
                x=self.x[ridx].value,
                pert_starts=np.asarray(pert_starts),
                pert_ends=np.asarray(pert_ends),
                ridx=ridx,
                subjname=subj.name,
                calculate_qpcr_loglik=calculate_qpcr_loglik,
                h5py_xname=self.x[ridx].name,
                target_acceptance_rate=self.target_acceptance_rate)
            if self.mp == 'debug':
                self.pool.append(worker)
            elif self.mp == 'full':
                pid = self.pool.add_worker(worker)
                self.worker_pids.append(pid)

        # Set the data to the latent values
        self.set_latent_as_data(update_values=False)

        self.total_n_datapoints = 0
        for ridx in range(self.G.data.n_replicates):
            self.total_n_datapoints += self.x[ridx].value.shape[0] * self.x[ridx].value.shape[1]

    def update(self):
        '''
        Copy-pasted from posterior.py implementation with some minor changes.
        '''
        self._strr = 'NA'
        if self.sample_iter < self.delay:
            return
        growth = self.G[STRNAMES.GROWTH_VALUE].value.ravel()
        self_interactions = self.G[STRNAMES.SELF_INTERACTION_VALUE].value.ravel()
        pv = self.G[STRNAMES.PROCESSVAR].value
        interactions = self.G[STRNAMES.INTERACTIONS_OBJ].get_datalevel_value_matrix(set_neg_indicators_to_nan=False)
        perts = None
        if self._there_are_perturbations:
            perts = []
            for perturbation in self.G.perturbations:
                perts.append(perturbation.item_array().reshape(-1,1))
            perts = np.hstack(perts)

        zero_inflation = [self.G[STRNAMES.ZERO_INFLATION].value[ridx] for ridx in range(self.G.data.n_replicates)]
        qpcr_vars = []
        for aaa in self.G[STRNAMES.QPCR_VARIANCES].value:
            qpcr_vars.append(aaa.value)

        kwargs = {'growth':growth, 'self_interactions':self_interactions,
            'pv':pv, 'interactions':interactions, 'perturbations':perts,
            'zero_inflation_data': zero_inflation, 'qpcr_variances':qpcr_vars}

        for ridx in range(self.G.data.n_replicates):
            _, x, acc_rate = self.pool[ridx].persistent_run(**kwargs)
            self.x[ridx].value = x

        self.set_latent_as_data()


class SubjectLogTrajectorySetWithoutQPCR(posterior.SubjectLogTrajectorySetMP):
    def __init__(
            self,
            biomass_mean: float,
            biomass_cv: float,
    ):
        """
        - "Biomass mean" is the desired mean for the prior, modeled as a lognormal distribution. calculate the appropriate normal mean and variance to reach the target lognormal mean/var.
        - "cv" = "coefficient of variation" = stdev / mean
        """
        self.biomass_lognormal_loc = np.log(biomass_mean)
        # note: the formula for CV for a lognormal RV X is: CV[X] = sqrt(exp(\sigma^2) - 1).
        self.biomass_lognormal_scale = np.sqrt(np.log(1 + np.square(biomass_cv)))
        super().__init__()

    def data_loglik_wo_intermediates(self) -> float:
        '''data loglikelihood with intermediate timepoints
        '''
        sum_q = self.sum_q[self.tidx]  # Sum of abundance latent state x_t
        log_sum_q = np.log(sum_q)
        rel = self.curr_x[self.tidx] / sum_q

        try:
            negbin = negbin_loglikelihood_MH_condensed(
                k=self.curr_reads,
                m=self.curr_read_depth * rel,
                dispersion=self.a0/rel + self.a1)
        except:
            negbin = negbin_loglikelihood_MH_condensed_not_fast(
                k=self.curr_reads,
                m=self.curr_read_depth * rel,
                dispersion=self.a0/rel + self.a1)

        biomass_ll = pl.random.normal.logpdf(value=log_sum_q, loc=self.biomass_lognormal_loc, scale=self.biomass_lognormal_scale)
        return negbin + biomass_ll


def initialize_graph_no_qpcr(
        params: config.MDSINE2ModelConfig,
        graph_name: str,
        subjset: Study,
        biomass_mean: float,
) -> BaseMCMC:
    '''
    Builds the probabilistic model for the MDSINE2 without time-series qPCR measurements.
    This is an adaptation of initialize_graph() from run.py in MDSINE2 implementation.
    :return: BaseMCMC object, which can then later invoke run_graph().
    '''
    # ----------
    # Type Check
    # ----------
    GRAPH = pl.Graph(name=graph_name, seed=params.SEED)
    GRAPH.as_default()
    pl.seed(params.SEED)

    # -----------------
    # Make the basepath
    # -----------------
    basepath = params.OUTPUT_BASEPATH
    os.makedirs(basepath, exist_ok=True)

    logger.info("Replacing all qPCR numbers with {} prior to renormalization".format(biomass_mean))
    for subj in subjset:
        for t in subj.qpcr:
            # Override qPCR measurements with constant values.
            qpcr_t = subj.qpcr[t]
            qpcr_t = qPCRdata(cfus=biomass_mean * np.ones(qpcr_t._raw_data.shape), mass=1., dilution_factor=1.)
            subj.qpcr[t] = qpcr_t

    # Normalize the qpcr measurements for numerical stability
    # -------------------------------------------------------
    if params.QPCR_NORMALIZATION_MAX_VALUE is not None:
        subjset.normalize_qpcr(max_value=params.QPCR_NORMALIZATION_MAX_VALUE)
        logger.info('Normalizing abundances for a max value of {}. Normalization ' \
                    'constant: {:.4E}'.format(params.QPCR_NORMALIZATION_MAX_VALUE,
                                              subjset.qpcr_normalization_factor))

        params.INITIALIZATION_KWARGS[STRNAMES.FILTERING]['v2'] *= subjset.qpcr_normalization_factor
        params.INITIALIZATION_KWARGS[STRNAMES.SELF_INTERACTION_VALUE]['rescale_value'] = \
            subjset.qpcr_normalization_factor
        biomass_mean = biomass_mean * subjset.qpcr_normalization_factor

    # ---------------------------------
    # Instantiate the posterior classes
    # ---------------------------------
    taxa = subjset.taxa
    data = design_matrices.Data(subjects=subjset, G=GRAPH,
                             zero_inflation_transition_policy=params.ZERO_INFLATION_TRANSITION_POLICY)
    clustering = Clustering(clusters=None, items=taxa, G=GRAPH, name=STRNAMES.CLUSTERING_OBJ)

    # Interactions
    var_interactions = posterior.PriorVarInteractions(
        prior=pl.variables.SICS(
            dof=pl.Constant(None, G=GRAPH),
            scale=pl.Constant(None, G=GRAPH),
            G=GRAPH), G=GRAPH)
    mean_interactions = posterior.PriorMeanInteractions(
        prior=pl.variables.Normal(
            loc=pl.Constant(None, G=GRAPH),
            scale2=pl.Constant(None, G=GRAPH),
            G=GRAPH), G=GRAPH)
    interaction_value = pl.variables.Normal(
        loc=mean_interactions, scale2=var_interactions, G=GRAPH)

    interaction_indicator = posterior.ClusterInteractionIndicatorProbability(
        prior=pl.variables.Beta(a=pl.Constant(None, G=GRAPH), b=pl.Constant(None, G=GRAPH), G=GRAPH),
        G=GRAPH)
    interactions = posterior.ClusterInteractionValue(
        prior=interaction_value, clustering=clustering, G=GRAPH)
    Z = posterior.ClusterInteractionIndicators(prior=interaction_indicator, G=GRAPH)

    # Growth
    var_growth = posterior.PriorVarMH(
        prior=pl.variables.SICS(
            dof=pl.Constant(None, G=GRAPH),
            scale=pl.Constant(None, G=GRAPH), G=GRAPH),
        child_name=STRNAMES.GROWTH_VALUE, G=GRAPH)
    mean_growth = posterior.PriorMeanMH(
        prior=pl.variables.TruncatedNormal(
            loc=pl.Constant(None, G=GRAPH),
            scale2=pl.Constant(None, G=GRAPH), G=GRAPH),
        child_name=STRNAMES.GROWTH_VALUE, G=GRAPH)
    prior_growth = pl.variables.TruncatedNormal(
        loc=mean_growth, scale2=var_growth,
        name='prior_{}'.format(STRNAMES.GROWTH_VALUE), G=GRAPH)
    growth = posterior.Growth(prior=prior_growth, G=GRAPH)

    # Self-Interactions
    var_si = posterior.PriorVarMH(
        prior=pl.variables.SICS(
            dof=pl.Constant(None, G=GRAPH),
            scale=pl.Constant(None, G=GRAPH), G=GRAPH),
        child_name=STRNAMES.SELF_INTERACTION_VALUE, G=GRAPH)
    mean_si = posterior.PriorMeanMH(
        prior=pl.variables.TruncatedNormal(
            loc=pl.Constant(None, G=GRAPH),
            scale2=pl.Constant(None, G=GRAPH), G=GRAPH),
        child_name=STRNAMES.SELF_INTERACTION_VALUE, G=GRAPH)
    prior_si = pl.variables.TruncatedNormal(
        loc=mean_si, scale2=var_si,
        name='prior_{}'.format(STRNAMES.SELF_INTERACTION_VALUE), G=GRAPH)
    self_interactions = posterior.SelfInteractions(prior=prior_si, G=GRAPH)

    # Process Variance
    prior_processvar = pl.variables.SICS(
        dof=pl.Constant(None, G=GRAPH),
        scale=pl.Constant(None, G=GRAPH), G=GRAPH)
    processvar = posterior.ProcessVarGlobal(G=GRAPH, prior=prior_processvar)

    # Clustering
    prior_concentration = pl.variables.Gamma(
        shape=pl.Constant(None, G=GRAPH),
        scale=pl.Constant(None, G=GRAPH),
        G=GRAPH)
    concentration = posterior.Concentration(
        prior=prior_concentration, G=GRAPH)
    cluster_assignments = posterior.ClusterAssignments(
        clustering=clustering, concentration=concentration,
        G=GRAPH, mp=params.MP_CLUSTERING)

    # Filtering and zero inflation
    filtering = FilteringWithoutQPCR(
        biomass_mean=biomass_mean,
        biomass_cv=0.10,  # make biomass standard deviation 10% of the mean.
        G=GRAPH,
        mp=params.MP_FILTERING,
        zero_inflation_transition_policy=params.ZERO_INFLATION_TRANSITION_POLICY
    )
    zero_inflation = posterior.ZeroInflation(zero_inflation_data_path=params.ZERO_INFLATION_DATA_PATH, G=GRAPH)

    # Perturbations
    if subjset.perturbations is not None:
        for pidx, subj_pert in enumerate(subjset.perturbations):
            name = subj_pert.name
            perturbation = ClusterPerturbationEffect(
                starts=subj_pert.starts, ends=subj_pert.ends,
                probability=pl.variables.Beta(
                    name=name + '_probability', G=GRAPH, value=None, a=None, b=None),
                clustering=clustering, G=GRAPH, name=name,
                signal_when_clusters_change=False, signal_when_item_assignment_changes=False)

            magnitude_var = posterior.PriorVarPerturbationSingle(
                prior=pl.variables.SICS(
                    dof=pl.Constant(None, G=GRAPH),
                    scale=pl.Constant(None, G=GRAPH), G=GRAPH),
                perturbation=perturbation, G=GRAPH)
            magnitude_mean = posterior.PriorMeanPerturbationSingle(
                prior=pl.variables.Normal(
                    loc=pl.Constant(None, G=GRAPH),
                    scale2=pl.Constant(None, G=GRAPH),
                    G=GRAPH),
                perturbation=perturbation, G=GRAPH)
            prior_magnitude = pl.variables.Normal(G=GRAPH, loc=magnitude_mean, scale2=magnitude_var)
            perturbation.magnitude.add_prior(prior_magnitude)

            prior_prob = pl.variables.Beta(
                a=pl.Constant(None, G=GRAPH),
                b=pl.Constant(None, G=GRAPH),
                G=GRAPH)
            perturbation.probability.add_prior(prior_prob)

        magnitude_var_perts = posterior.PriorVarPerturbations(G=GRAPH)
        magnitude_mean_perts = posterior.PriorMeanPerturbations(G=GRAPH)
        magnitude_perts = posterior.PerturbationMagnitudes(G=GRAPH)
        indicator_perts = posterior.PerturbationIndicators(G=GRAPH, need_to_trace=False, relative=True)
        indicator_prob_perts = posterior.PerturbationProbabilities(G=GRAPH)
    else:
        magnitude_perts = None
        pert_ind = None
        pert_ind_prob = None

    beta = posterior.GLVParameters(
        growth=growth, self_interactions=self_interactions,
        interactions=interactions, pert_mag=magnitude_perts, G=GRAPH)

    # Set qPCR variance priors and hyper priors
    qpcr_variances = posterior.qPCRVariances(G=GRAPH, L=params.N_QPCR_BUCKETS)

    # Set up inference and the inference order.
    # -----------------------------------------
    mcmc = pl.BaseMCMC(burnin=params.BURNIN, n_samples=params.N_SAMPLES, graph=GRAPH)
    order = []
    for name in params.INFERENCE_ORDER:
        if name == STRNAMES.QPCR_DOFS:
            continue
        if name == STRNAMES.QPCR_SCALES:
            continue
        if name == STRNAMES.QPCR_VARIANCES:
            continue
        if params.LEARN[name]:
            if not STRNAMES.is_perturbation_param(name):
                order.append(name)
            elif subjset.perturbations is not None:
                order.append(name)
    mcmc.set_inference_order(order)

    # Initialize the posterior and instantiate the design matrices
    # ------------------------------------------------------------
    for name in params.INITIALIZATION_ORDER:
        logger.info('Initializing {}'.format(name))
        if STRNAMES.is_perturbation_param(name) and subjset.perturbations is None:
            logger.info('Skipping over {} because it is a perturbation parameter ' \
                        'and there are no perturbations'.format(name))
            continue

        if name not in GRAPH:
            logger.warn("Node `{}` not found in graph, but was specified in parameters. Skipping.".format(name))
            continue

        # Call `initialize`
        try:
            GRAPH[name].initialize(**params.INITIALIZATION_KWARGS[name])
        except Exception as error:
            logger.critical('Initialization in `{}` failed with the parameters: {}'.format(
                name, params.INITIALIZATION_KWARGS[name]) + ' with the follwing error:\n{}'.format(
                error))
            for a in GRAPH._persistent_pntr:
                a.kill()
            raise

        # Initialize data matrices if necessary
        if name == STRNAMES.ZERO_INFLATION:
            # Initialize the basic data matrices after initializing filtering
            lhs = design_matrices.LHSVector(G=GRAPH, name='lhs_vector')
            lhs.build()
            growthDM = design_matrices.GrowthDesignMatrix(G=GRAPH)
            growthDM.build_without_perturbations()
            selfinteractionsDM = design_matrices.SelfInteractionDesignMatrix(G=GRAPH)
            selfinteractionsDM.build()
        if name == STRNAMES.CLUSTER_INTERACTION_INDICATOR:
            # Initialize the interactions data matrices after initializing the interactions
            interactionsDM = design_matrices.InteractionsDesignMatrix(G=GRAPH)
            interactionsDM.build()
        if name == STRNAMES.PERT_INDICATOR and subjset.perturbations is not None:
            # Initialize the perturbation data matrices after initializing the perturbations
            perturbationsDM = design_matrices.PerturbationDesignMatrix(G=GRAPH)
            perturbationsDM.base.build()
            perturbationsDM.M.build()
        if name == STRNAMES.PERT_VALUE and subjset.perturbations is not None:
            data.design_matrices[STRNAMES.GROWTH_VALUE].build_with_perturbations()

    # Setup filenames
    # ---------------
    hdf5_filename = os.path.join(basepath, config.HDF5_FILENAME)
    mcmc_filename = os.path.join(basepath, config.MCMC_FILENAME)
    mcmc.set_tracer(filename=hdf5_filename, checkpoint=params.CHECKPOINT)
    mcmc.set_save_location(mcmc_filename)
    params.save(os.path.join(basepath, config.PARAMS_FILENAME))
    pl.seed(params.SEED)

    return mcmc
