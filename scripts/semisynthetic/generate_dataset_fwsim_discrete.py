from collections import namedtuple
from pathlib import Path
from typing import Tuple, Dict, List, Set
import pickle as pkl

import argparse

import numba
import numpy as np
import scipy
import scipy.stats

import mdsine2 as md2
from mdsine2.names import STRNAMES
import matplotlib.pyplot as plt


# ======= Helper functions
def create_synthetic_dataset(
        source_study: md2.Study,
        forward_sims: Dict[str, np.ndarray],
        sim_read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator,
) -> md2.Study:
    print("Sampling.")
    read_samples_per_subj, qpcr_samples_per_subj = sample_data_from_fwsim(
        source_study,
        forward_sims,
        sim_read_depth,
        negbin_a0,
        negbin_a1,
        qpcr_noise_scale,
        rng,
    )

    # create the study object.
    synthetic_study = md2.Study(taxa=source_study.taxa, name='simulated')
    for source_subj in source_study:
        synthetic_study.add_subject(name=source_subj.name)
        new_subj = synthetic_study[source_subj.name]
        for tidx, t in enumerate(source_subj.times):
            sampled_reads = read_samples_per_subj[source_subj.name]  # (n_taxa, n_timepoints)
            sampled_qpcr = qpcr_samples_per_subj[source_subj.name]
            new_subj.reads[t] = sampled_reads[:, tidx]
            new_subj.qpcr[t] = md2.qPCRdata(cfus=sampled_qpcr[tidx], mass=1., dilution_factor=1.)
        new_subj.times = source_subj.times
    synthetic_study.perturbations = source_study.perturbations
    print("synthetic study name: {}".format(synthetic_study.name))
    return synthetic_study


def create_synthetic_replicates(
        source_study: md2.Study,
        forward_sims: Dict[str, np.ndarray],
        source_subj_name: str,
        source_subj_timepoints: List[float],  # timepoints
        num_physical_replicates: int,
        replicate_study_name: str,
        negbin_a0: float,
        negbin_a1: float,
        num_reads: int,
        qpcr_noise_scale: float,
        rng: np.random.Generator
) -> md2.Study:
    # Extract the source subj and related data
    src_subj = source_study[source_subj_name]
    timepoint_indices = {t: t_idx for t_idx, t in enumerate(src_subj.times)}

    # Extract the "true" biomass (the posterior median of the filtered state abundances)
    true_traj = forward_sims[source_subj_name]

    # Make the study object
    replicate_study = md2.Study(taxa=source_study.taxa, name=replicate_study_name)
    for t in source_subj_timepoints:
        replicate_study.add_subject(name=f'M{source_subj_name}-D{t}')

    # Add times for each subject
    replicate_times = np.arange(0., num_physical_replicates, 1.0)
    for replicate_subj in replicate_study:
        replicate_subj.times = replicate_times

    for replicate_subj, subj_t in zip(replicate_study, source_subj_timepoints):
        # Extract timepoint abundances
        if subj_t not in timepoint_indices:
            raise Exception("Couldn't find timepoint {} in subject {}.".format(subj_t, source_subj_name))
        t_idx = timepoint_indices[subj_t]

        conc = true_traj[:, t_idx]
        total_mass = np.sum(conc)
        rel_mass = conc / total_mass

        # simulate physical replicates
        for repl_t in replicate_times:
            # Make the reads
            replicate_subj.reads[repl_t] = np.array([
                negative_binomial(phi=ra_ratio * num_reads, eps=negbin_a1 + (negbin_a0 / ra_ratio), rng=rng)
                if ra_ratio > 0
                else 0
                for ra_ratio in rel_mass
            ])

            # simulate technical/measurement triplicates
            triplicates = np.exp(
                np.log(total_mass) +
                rng.normal(loc=0.0, scale=qpcr_noise_scale, size=3)
            )
            replicate_subj.qpcr[repl_t] = md2.qPCRdata(cfus=triplicates, mass=1., dilution_factor=1.)
    print("replicate study name: {}".format(replicate_study.name))
    return replicate_study



GLVParamSet = namedtuple(
    'GLVParamSet',
    [
        'growth', 'interactions', 'perturbations'
    ]
)


def extract_glv_model(
        mcmc: md2.BaseMCMC,
        study: md2.Study,
) -> Tuple[GLVParamSet, np.ndarray, Dict[str, np.ndarray]]:
    mcmc_log_likelihoods = evaluate_likelihoods(mcmc, study)
    best_idx = np.argmax(mcmc_log_likelihoods)


    glv_params = GLVParamSet(
        growth=growths[i],
        interactions=interactions[i],
        perturbations=[p[i] for p in perturbations]
    )
    return glv_params, coclusterings, fwsim


def clustering_sizes(cluster_assignments) -> np.ndarray:
    cluster_ids = set(cluster_assignments)
    return np.array([
        np.sum(cluster_ids == c_id)
        for c_id in cluster_ids
    ])


# ======================================== LIKELIHOOD FUNCTIONS/HELPERS
def crp_log_likelihood(alpha, clustering):
    sizes = clustering_sizes(clustering)
    ll = 0.0
    cumulative_sz = 0
    for n in sizes:
        # numerator: alpha * 1 * 2 * ... * (n-1)
        ll += np.sum(np.log(alpha) + np.log(np.arange(1, n, step=1)))

        # denominator: (C+alpha) + (C+alpha+1) + ... + (C+n-1+alpha)
        ll -= np.sum(np.log(cumulative_sz + alpha + np.arange(0, n)))
    return ll


def module_interaction_log_likelihood(clustering: np.ndarray, interactions: np.ndarray, indicator_prior_a, indicator_prior_b, strength_prior_means, strength_prior_scales) -> float:
    clustering_ids = set(clustering)
    c_rep_indices = [
        np.argmax(clustering == c_id)
        for c_id in clustering_ids
    ]
    submatrix = interactions[c_rep_indices, :][:, c_rep_indices]

    # First, evaluate the zero-locations (log-likelihood of the indicator being zero)
    n_zeros = np.sum(np.isnan(submatrix))
    zero_ll = n_zeros * scipy.stats.betabinom.logpmf(n=1, k=0, a=indicator_prior_a, b=indicator_prior_b)

    # Next, evaluate the nonzero-locations (LL of indicator being nonzero, plus LL of strength)
    nonzeros = submatrix[~np.isnan(submatrix)]
    nonzero_ll = nonzeros.size * scipy.stats.betabinom.logpmf(n=1, k=1, a=indicator_prior_a, b=indicator_prior_b)
    strength_ll = np.sum(scipy.stats.norm.logpdf(loc=strength_prior_means, scale=strength_prior_scales, x=nonzeros))
    return zero_ll + nonzero_ll + strength_ll


def pert_log_likelihood(pert_str, clustering, pert_indicator_a, pert_indicator_b, pert_strength_mean, pert_strength_scale):
    clustering_ids = set(clustering)
    c_rep_indices = [
        np.argmax(clustering == c_id)
        for c_id in clustering_ids
    ]
    subsection = pert_str[c_rep_indices]

    # First, evaluate the zero-locations (LL of indicator being zero)
    n_zeros = np.sum(np.isnan(subsection))
    zero_ll = n_zeros * scipy.stats.betabinom.logpmf(n=1, k=0, a=pert_indicator_a, b=pert_indicator_b)

    # Next, evaluate the nonzero-locations (LL of indicator being nonzero, plus LL of strength)
    nonzeros = subsection[~np.isnan(subsection)]
    nonzero_ll = nonzeros.size * scipy.stats.betabinom.logpmf(n=1, k=1, a=pert_indicator_a, b=pert_indicator_b)
    strength_ll = np.sum(scipy.stats.norm.logpdf(loc=pert_strength_mean, scale=pert_strength_scale, x=nonzeros))
    return zero_ll + nonzero_ll + strength_ll


def forward_process_ll(
        subj: md2.Subject,
        subj_latent_x: np.ndarray,
        growths: np.ndarray,
        interactions: np.ndarray,
        all_pert_strs: Dict[str, np.ndarray],
        process_vars: np.ndarray
) -> np.ndarray:
    # This is vectorized.
    n_samples = subj_latent_x.shape[0]
    assert growths.shape[0] == n_samples
    assert interactions.shape[0] == n_samples
    for p in all_pert_strs.values():
        assert p.shape[0] == n_samples
    assert process_vars.shape[0] == n_samples

    # Formulate the perturbation on/off calculation per timepoint.
    pert_order = sorted(all_pert_strs.keys())
    perturbation_signals = np.stack([
        (subj.times[:-1] >= subj.parent.perturbations[pert_name].starts[subj.name])
        &
        (subj.times[:-1] <= subj.parent.perturbations[pert_name].ends[subj.name])
        for pert_name in pert_order
    ], axis=0)  # shape is (N_perts, N_times-1), indicates whether perturbation is applied at a particular timepoint.
    pert_strs_matrix = np.stack([
        all_pert_strs[pert_name]
        for pert_name in pert_order
    ], axis=2)  # shape is (N_samples, N_taxa, N_perts)

    # latent process likelihood
    dt = np.diff(subj.times)
    subj_latent_mu_log = (
            np.log(subj_latent_x[:, :, :-1])  # shape is (N_samples, N_taxa, N_times-1)
            + (
                    ( # r+Ax is implemented here
                            ( # r is implemented here: r = (growth)*(1+pert)
                                    growths  # shape is (N_samples, N_taxa)
                                    * (1 + pert_strs_matrix @ perturbation_signals)  # shape is (N_samples, N_taxa, N_perts) @ (N_perts, N_times-1) --> (N_samples, N_taxa, N_times-1)
                            )
                            +
                            ( # A is implemented here: A_ii = self-interaction, A_ij = interaction.
                                # Note: The "@" operator auto-broadcasts (see np.matmul).
                                    interactions  # shape is (N_samples, N_taxa, N_taxa)
                                    @ subj_latent_x[:, :, :-1]  # shape is (N_samples, N_taxa, N_times-1)
                            )
                    )
                    * dt[None, None, :]  # shape is (1, 1, N_times-1)
            )
    )  # final shape is (N_samples, N_taxa, N_times-1)

    # raw logpdf shape is (N_samples, N_taxa, N_times-1), so invoke sum twice at the end to reduce to (N_samples)
    return scipy.stats.norm.logpdf(loc=subj_latent_mu_log, scale=process_vars[:, None, None], x=subj_latent_x[:, :, 1:]).sum(axis=-1).sum(axis=-1)


def filtered_data_ll(
        subj: md2.Subject,
        subj_latent_x: np.ndarray,
        read_negbin_d0: float,
        read_negbin_d1: float
) -> np.ndarray:
    reads_matrix = np.stack([subj.reads[t] for t in subj.times], axis=-1)  # shape is (N_taxa, N_times)
    n_taxa = reads_matrix.shape[0]
    n_times = len(subj.times)

    assert subj_latent_x.shape[1] == n_taxa
    assert subj_latent_x.shape[2] == n_times

    # ===== READ LLS
    total_reads = reads_matrix.sum(axis=0)  # shape is (N_times)

    rel_abunds = subj_latent_x / np.sum(subj_latent_x, axis=1, keepdims=True)  # (N_samples, N_taxa, N_times)
    nb_phi = total_reads[None, None, :] * rel_abunds
    nb_eps = read_negbin_d1 + read_negbin_d0 / rel_abunds

    nb_p = np.reciprocal(1 + nb_phi * nb_eps)
    nb_n = np.reciprocal(nb_eps)
    reads_ll = scipy.stats.nbinom.logpmf(n=nb_n, p=nb_p, k=reads_matrix).sum(axis=-1).sum(axis=-1)  # LL shape is (N_samples, N_taxa, N_times), then summed into (N_samples)

    # ===== qPCR LLS
    log_qpcr_matrix = np.stack([
        np.log(subj.qpcr[t]._raw_data)
        for t in subj.times
    ], axis=0)  # shape is (N_times, N_replicates)
    qpcr_mu = np.log(np.sum(subj_latent_x, axis=1))  # shape is (N_samples, N_times)
    qpcr_scale = np.std(log_qpcr_matrix, axis=1)  # shape is (N_times)

    # target shape is (N_samples, N_times, N_replicates)
    qpcr_ll = scipy.stats.norm.logpdf(
        loc=np.expand_dims(qpcr_mu, axis=2),
        scale=qpcr_scale[None, :, None],
        x=np.expand_dims(log_qpcr_matrix, axis=0)
    ).sum(axis=-1).sum(axis=-1)

    # ===== Output
    return reads_ll + qpcr_ll


def evaluate_likelihoods(
        mcmc: md2.BaseMCMC, study: md2.Study
):
    # Growth prior ll
    growths = mcmc.graph[STRNAMES.GROWTH].get_trace_from_disk(section='posterior')
    growth_strength_means = mcmc.graph[STRNAMES.PRIOR_MEAN_GROWTH].get_trace_from_disk(section='posterior')
    growth_strength_vars = mcmc.graph[STRNAMES.PRIOR_VAR_GROWTH].get_trace_from_disk(section='posterior')
    growth_strength_prior_dof = mcmc.graph[STRNAMES.PRIOR_VAR_GROWTH].prior.dof.value
    growth_strength_prior_tau2 = mcmc.graph[STRNAMES.PRIOR_VAR_GROWTH].prior.scale.value
    growth_prior_lls = (
        scipy.stats.norm.logpdf(loc=mcmc.graph[STRNAMES.PRIOR_MEAN_GROWTH].prior.loc.value, scale=mcmc.graph[STRNAMES.PRIOR_MEAN_GROWTH].prior.scale, x=growth_strength_means)
        + scipy.stats.invgamma.logpdf(a=(0.5 * growth_strength_prior_dof), scale=(0.5 * growth_strength_prior_dof * growth_strength_prior_tau2), x=growth_strength_vars)
    )
    growth_lls = scipy.stats.norm.logpdf(loc=growth_strength_means, scale=np.sqrt(growth_strength_vars), x=growths)

    # Self interaction prior ll
    self_interactions = mcmc.graph[STRNAMES.SELF_INTERACTION_VALUE].get_trace_from_disk(section='posterior')
    self_interaction_strength_means = mcmc.graph[STRNAMES.PRIOR_MEAN_SELF_INTERACTIONS].get_trace_from_disk(section='posterior')
    self_interaction_strength_vars = mcmc.graph[STRNAMES.PRIOR_VAR_SELF_INTERACTIONS].get_trace_from_disk(section='posterior')
    self_interaction_strength_prior_dof = mcmc.graph[STRNAMES.PRIOR_VAR_SELF_INTERACTIONS].prior.dof.value
    self_interaction_strength_prior_tau2 = mcmc.graph[STRNAMES.PRIOR_VAR_SELF_INTERACTIONS].prior.scale.value
    self_interaction_prior_lls = (
        scipy.stats.norm.logpdf(loc=mcmc.graph[STRNAMES.PRIOR_MEAN_SELF_INTERACTIONS].prior.loc.value, scale=mcmc.graph[STRNAMES.PRIOR_MEAN_SELF_INTERACTIONS].prior.scale, x=self_interaction_strength_means)
        + scipy.stats.invgamma.logpdf(a=(0.5 * self_interaction_strength_prior_dof), scale=(0.5 * self_interaction_strength_prior_dof * self_interaction_strength_prior_tau2), x=self_interaction_strength_vars)
    )
    self_interaction_lls = scipy.stats.norm.logpdf(loc=self_interaction_strength_means, scale=np.sqrt(self_interaction_strength_vars), x=self_interactions)

    # clustering prior ll
    alphas = mcmc.graph[STRNAMES.CONCENTRATION].get_trace_from_disk(section='posterior')
    clusterings = np.stack(list(extract_clusterings(mcmc)))

    clustering_lls = scipy.stats.gamma.logpdf(a=1e5, scale=1e-5, x=alphas)
    clustering_lls += np.array([
        crp_log_likelihood(alpha, clustering) for alpha, clustering in zip(alphas, clusterings)
    ])

    # interaction prior ll
    interactions = mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk(section='posterior')
    interaction_indicator_a = mcmc.graph[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB].prior.a.value
    interaction_indicator_b = mcmc.graph[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB].prior.b.value
    interaction_strength_means = mcmc.graph[STRNAMES.PRIOR_MEAN_INTERACTIONS].get_trace_from_disk(section='posterior')
    interaction_strength_variances = mcmc.graph[STRNAMES.PRIOR_VAR_INTERACTIONS].get_trace_from_disk(section='posterior')
    interaction_strength_prior_dof = mcmc.graph[STRNAMES.PRIOR_VAR_INTERACTIONS].prior.dof.value
    interaction_strength_prior_tau2 = mcmc.graph[STRNAMES.PRIOR_VAR_INTERACTIONS].prior.scale.value
    interaction_prior_lls = (
        scipy.stats.norm.logpdf(loc=mcmc.graph[STRNAMES.PRIOR_MEAN_INTERACTIONS].prior.loc.value, scale=mcmc.graph[STRNAMES.PRIOR_VAR_INTERACTIONS].prior.scale, x=interaction_strength_means)
        + scipy.stats.invgamma.logpdf(a=(0.5 * interaction_strength_prior_dof), scale=(0.5 * interaction_strength_prior_dof * interaction_strength_prior_tau2), x=interaction_strength_variances)
    )
    interaction_lls = np.array([
        module_interaction_log_likelihood(
            clustering, interaction_matrix,
            interaction_indicator_a,
            interaction_indicator_b,
            interaction_strength_means,
            np.sqrt(interaction_strength_variances)
        )
        for clustering, interaction_matrix in zip(clusterings, interactions)
    ])

    # Perturbation prior ll
    pert_indicator_a = 0.5   # weak-agnostic
    pert_indicator_b = 0.5   # weak-agnostic

    all_pert_strs = {}
    perts_prior_lls = np.zeros(growth_lls.shape, dtype=float)
    perts_ll = np.zeros(growth_lls.shape, dtype=float)

    for pert in mcmc.graph.perturbations:
        pert_strengths = pert.magnitude.get_trace_from_disk(section='posterior')
        all_pert_strs[pert.name] = pert_strengths

        pert_strength_means = pert.magnitude.prior.loc.get_trace_from_disk(section='posterior')
        pert_strength_variances = pert.magnitude.prior.scale2.get_trace_from_disk(section='posterior')
        perts_prior_lls += scipy.stats.norm.logpdf(loc=pert.magnitude.prior.loc.value, scale=pert.magnitude.prior.scale, x=pert_strength_means)
        perts_prior_dof = pert.magnitude.prior.dof.value
        perts_prior_tau2 = pert.magnitude.prior.scale.value
        perts_prior_lls += scipy.stats.invgamma.logpdf(a=(0.5 * perts_prior_dof), scale=(0.5 * perts_prior_dof * perts_prior_tau2), x=pert_strength_variances)
        perts_ll += np.array([
            pert_log_likelihood(pert_str, clustering, pert_indicator_a, pert_indicator_b, pert_mean, pert_scale)
            for pert_str, clustering, pert_mean, pert_scale in zip(
                pert_strengths,
                clusterings,
                pert_strength_means,
                np.sqrt(pert_strength_variances)
            )
        ])

    # Data ll, given parameters
    interactions[np.isnan(interactions)] = 0.
    for i in range(self_interactions.shape[1]):
        interactions[:, i, i] = -np.absolute(self_interactions[:, i])
    for pert_name, pert_str in all_pert_strs.items():
        pert_str[np.isnan(pert_str)] = 0.

    data_ll = np.zeros(growth_lls.shape, dtype=float)
    n_samples = growths.shape[0]
    n_taxa = len(study.taxa)
    process_vars = mcmc.graph[STRNAMES.PROCESSVAR].get_trace_from_disk(section='posterior')
    read_negbin_d0 = mcmc.graph[STRNAMES.FILTERING].a0
    read_negbin_d1 = mcmc.graph[STRNAMES.FILTERING].a1
    for subj, subj_traj_variable in zip(study, mcmc.graph[STRNAMES.FILTERING].x.value):
        n_timepoints = len(subj.times)
        subj_latent_x = subj_traj_variable.get_trace_from_disk(section='posterior')  # shape is (samples, n_taxa, n_timepoints)
        assert subj_latent_x.shape[0] == n_samples
        assert subj_latent_x.shape[1] == n_taxa
        assert subj_latent_x.shape[2] == n_timepoints
        subj_filter_ll = forward_process_ll(subj, subj_latent_x, growths, interactions, all_pert_strs, process_vars)
        subj_data_ll = filtered_data_ll(subj, subj_latent_x, read_negbin_d0, read_negbin_d1)

        assert subj_filter_ll.size == n_samples
        assert subj_data_ll.size == n_samples

        data_ll += subj_filter_ll + subj_data_ll
    return (
            growth_prior_lls + growth_lls
            + self_interaction_prior_lls + self_interaction_lls
            + clustering_lls
            + interaction_prior_lls + interaction_lls
            + perts_prior_lls + perts_ll
            + data_ll
    )
# ======================================== END: LIKELIHOOD FUNCTIONS/HELPERS


def forward_simulate(
        glv_params: GLVParamSet,
        study: md2.Study,
        sim_max: float
) -> Dict[str, np.ndarray]:
    """Forward simulation for all subjects in a Study object"""
    return {
        subj.name: forward_simulate_subject(glv_params, study, subj, sim_max)[0]
        for subj in study
    }


def forward_simulate_subject(
        glv_params: GLVParamSet,
        study: md2.Study,
        subject: md2.Subject,
        sim_dt: float,
        sim_max: float,
        time_subset: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """Forward simulation for a single subject"""
    # ======= Perturbations
    perturbations_start = []
    perturbations_end = []
    if study.perturbations is not None:
        for pert in study.perturbations:
            perturbations_start.append(pert.starts[subject.name])
            perturbations_end.append(pert.ends[subject.name])

    dyn = md2.model.gLVDynamicsSingleClustering(
        growth=glv_params.growth,
        interactions=glv_params.interactions,
        perturbations=glv_params.perturbations,
        perturbation_starts=perturbations_start,
        perturbation_ends=perturbations_end,
        start_day=subject.times[0],
        sim_max=sim_max
    )

    # print("Data shape: ", subject.matrix()['abs'].shape)
    initial_conditions = subject.matrix()['abs'][:, 0] + 1e4
    if time_subset:
        x = md2.integrate(
            dynamics=dyn,
            initial_conditions=np.expand_dims(initial_conditions, 1),
            dt=sim_dt,
            final_day=subject.times[-1],
            subsample=True,
            times=subject.times
        )
    else:
        x = md2.integrate(
            dynamics=dyn,
            initial_conditions=np.expand_dims(initial_conditions, 1),
            dt=sim_dt,
            final_day=subject.times[-1],
            subsample=False
        )
    fwsim_values = x['X']
    times = x['times']
    return fwsim_values, times


# def integrate_fwsim_discrete_euler(
#         growths: np.ndarray,
#         interactions: np.ndarray,
#         perts,
#         pert_starts,
#         pert_ends,
#         x0,
#         times,
#         pert_is_multplicative: bool,
#         sim_max: float,
# ) -> np.ndarray:
#     @numba.njit
#     def grad(log_x: np.ndarray, g: np.ndarray, A: np.ndarray):
#         return g + A @ np.exp(log_x)
#
#     sim_max_log = np.log(sim_max)
#     n_taxa = len(growths)
#     n_times = len(times)
#     log_x = np.zeros(shape=(n_taxa, n_times))
#     log_x[:, 0] = np.log(x0)
#
#     for t_idx in range(n_times - 1):
#         t = times[t_idx]
#         dt = times[t_idx + 1] - times[t_idx]
#         g = growths
#         for pert_strengths, pert_start, pert_end in zip(perts, pert_starts, pert_ends):
#             if t >= pert_start and t <= pert_end:
#                 if pert_is_multplicative:
#                     g = g * (1 + pert_strengths)
#                 else:
#                     g = g + pert_strengths
#         next_value = log_x[:, t_idx] + grad(log_x[:, t_idx], g, interactions) * dt
#         next_value[next_value > sim_max_log] = sim_max_log  # cap maximum value for numerical stability
#         log_x[:, t_idx + 1] = next_value
#     return np.exp(log_x)


# =========== helpers
# def extract_clustering_from_matrix(_mat) -> np.ndarray:
#     n_items = _mat.shape[0]
#     items_left = set(range(n_items))
#
#     clustering_assignments = np.zeros(n_items, dtype=int)
#
#     c_idx = -1
#     while len(items_left) > 0:
#         c_idx += 1  # new cluster
#         x = items_left.pop()  # member
#         clustering_assignments[x] = c_idx
#
#         # get all guys in the same cluster
#         to_remove = set()
#         for y in items_left:
#             if _mat[x,y]:
#                 clustering_assignments[y] = c_idx
#                 to_remove.add(y)
#         items_left = items_left.difference(to_remove)
#     return clustering_assignments


# def evaluate_parameter_fwsim(
#         params: GLVParamSet,
#         study: md2.Study,
#         sim_max: float,
#         coclust_matrix: np.ndarray
# ) -> Tuple[float, Set[int], Set[int], int]:
#     """
#     Evaluate the fwsim error for each subject (And sum them).
#     """
#     fwsims = forward_simulate(params, study, sim_max)
#     errors = []
#     eps = 1e-5
#     for subj_name, fwsim_trajs in fwsims.items():
#         subj = study[subj_name]
#         measurements = subj.matrix()['abs']
#
#         # RMS-log10 error, excluding the points where data says zero.
#         # fwsim_trajs: [taxa x timepoints] array
#         # measurements: [taxa x timepoints] array, with zeroes in some of the entries
#
#         mask = measurements > 0
#         diff_logs = np.log10(fwsim_trajs + eps) - np.log10(measurements + eps)
#         subj_err = np.sqrt(np.mean(np.square(diff_logs[mask])))
#         errors.append(subj_err)
#
#     # What taxa survived? (Only keep taxa which surpasses 1e5 abundance for some timepoint in some synthetic mouse)
#     taxa_max_abundance = np.stack([
#         fwsim_trajs.max(axis=-1)  # length=n_taxa, after maximizing across timepoints
#         for subj_name, fwsim_trajs in fwsims.items()
#     ], axis=0).max(axis=0)  # length=n_taxa, max across subjects
#
#     clust_assignments = extract_clustering_from_matrix(coclust_matrix)
#     leftover_module_set = set()
#     leftover_taxa_set = set()
#     for taxa_idx, taxa in enumerate(study.taxa):
#         if taxa_max_abundance[taxa_idx] > 1e5:
#             leftover_module_set.add(clust_assignments[taxa_idx])
#             leftover_taxa_set.add(taxa_idx)
#     total_modules = np.max(clust_assignments) + 1
#
#     print("leftover taxa: {}".format(len(leftover_taxa_set)))
#     return np.sum(errors), leftover_taxa_set, leftover_module_set, total_modules


def sample_data_from_fwsim(
        study: md2.Study,
        forward_sims: Dict[str, np.ndarray],
        read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        rng: np.random.Generator,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    :return: A tuple.
    (1) A dictionary of (Taxa x Timepoint) abundance samples per subject,
    (2) A dictionary of (Timepoint x Replicate) qPCR samples per subject.
    """
    reads_per_subject = {}
    qpcr_per_subject = {}
    for subj_idx, subj in enumerate(study):
        trajs = forward_sims[subj.name]

        if trajs.shape[0] != len(study.taxa):
            raise ValueError("Matrix's n_taxa dimension ({}) doesn't match taxaset size ({})".format(
                trajs.shape[0], len(study.taxa)
            ))
        if trajs.shape[1] != len(subj.times):
            raise ValueError("Matrix's n_times dimension ({}) doesn't match subj times size ({})".format(
                trajs.shape[1], len(subj.times)
            ))

        read_counts = sample_reads(read_depth, trajs, negbin_a0, negbin_a1, rng=rng)  # n_taxa x n_timepoints
        reads_per_subject[subj.name] = read_counts

        qpcrs = sample_qpcr(np.sum(trajs, axis=0), qpcr_noise_scale, n_replicates=3, rng=rng)
        qpcr_per_subject[subj.name] = qpcrs
    return reads_per_subject, qpcr_per_subject


def render_plots(study: md2.Study, sims: Dict[str, np.ndarray], plot_dir: Path):
    width = 6
    height_per_subj = 6
    fmt = 'pdf'
    for taxa_idx, taxa in enumerate(study.taxa):
        plot_path = plot_dir / f'{taxa.name}.{fmt}'

        fig, axes = plt.subplots(len(study), 1, figsize=(width, height_per_subj * len(study)))
        for subj, ax in zip(study, axes):
            measurements = subj.matrix()['abs'][taxa_idx, :]
            trajs = sims[subj.name]
            ax.plot(subj.times, trajs[taxa_idx, :])
            ax.plot(subj.times, measurements, marker='x', color='black', linestyle=':')
            ax.set_title(f'Subject {subj.name}')
            ax.set_yscale('log')
        fig.suptitle(f"Taxa {taxa.name}")
        plt.savefig(plot_path, format=fmt)
        plt.close(fig)


def sample_reads(
        total_reads: int,
        abundance_trajectories: np.ndarray,
        negbin_a0: float,
        negbin_a1: float,
        rng: np.random.Generator
) -> np.ndarray:
    """
    Sample reads from the (Taxa, Timepoints)-shaped array of trajectories.
    :param total_reads: The read depth to simulate.
    :param abundance_trajectories: The (Taxa, Timepoints)-shaped array of trajectories.
    :param dirichlet_alpha: The dirichlet "alpha" scaling parameter, defined as \alpha := \sum_i \alpha_i.
    :param rng: The numpy Generator object to use (for reproducibility with seeds)
    :return: Abundances in the shape (Taxa, Timepoints).
    """
    _, n_timepoints = abundance_trajectories.shape
    total_mass = np.sum(abundance_trajectories, axis=0)
    reads_all = []  # A list of (Taxa,)-shaped arrays.
    for tidx in range(n_timepoints):
        rel_abunds = abundance_trajectories[:, tidx] / total_mass[tidx]  # indexed by taxa

        phis = rel_abunds * total_reads
        epsilons = negbin_a1 + (negbin_a0 / rel_abunds)
        reads_all.append(
            np.array([
                negative_binomial(phi=phi, eps=eps, rng=rng)
                if ~np.isinf(eps)
                else 0
                for phi, eps in zip(phis, epsilons)
            ])
        )
    return np.stack(reads_all, axis=1)


def negative_binomial(phi, eps, rng: np.random.Generator) -> int:
    """
    :param phi: the mean parameter
    :param eps: the dispersion parameter
    :return:

    Refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html
    which gives:
     p = (mean) / (var) = 1 / (1 + phi * eps)
     n = (mean^2) / (var - mean) = 1 / eps
    """
    p = 1 / (1 + phi * eps)
    n = 1 / eps
    return rng.negative_binomial(n=n, p=p)


def sample_qpcr(
        total_mass: np.ndarray,
        noise_scale: float,
        n_replicates: int,
        rng: np.random.Generator
) -> np.ndarray:
    """
    Sample qPCR measurements for the specified biomass trajectory.
    :param total_mass: A 1-D array specifying biomass for each timepoint.
    :param noise_scale: The log-normal noise scale (std-deviation) parameter.
    :param n_replicates: The number of replicates (e.g. 3 for triplicates) of qpcr measurements.
    :param rng: The numpy Generator object to use (for reproducibility with seeds)
    """
    n_times, = total_mass.shape
    return np.exp(
        np.expand_dims(np.log(total_mass), axis=1)  # shape (T, 1)
        +
        noise_scale * rng.normal(loc=0.0, scale=1.0, size=(n_times, n_replicates))
    )


# ============ CLI interface
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--study', '-s', dest="study_path",
        type=str, required=True,
        help="The path to the original (real) datas Study pickle file."
    )
    parser.add_argument(
        '--growths', dest='growth_path',
        type=str, required=True, help="The path to the .npy growth rate MCMC samples."
    )
    parser.add_argument(
        '--interactions', dest='interactions_path',
        type=str, required=True, help="The path to the .npy interactions MCMC samples."
    )
    parser.add_argument(
        '--perts', dest='perts_path',
        type=str, required=True, help="The path to the .npz perturbation strengths MCMC samples."
    )
    parser.add_argument(
        '--coclust', dest='coclust_path',
        type=str, required=True, help="The path to the .npz coclustering MCMC samples."
    )
    parser.add_argument(
        '--truth-dir', '-t', dest='ground_truth_dir',
        type=str, required=True,
        help="A directory to store ground truth values."
    )
    parser.add_argument(
        '--out', '-o', dest='out_path',
        type=str, required=True,
        help="The desired full path to which the synthetic study object will be saved."
    )
    parser.add_argument(
        '--replicate-out', '-ro', dest='replicate_out_path',
        type=str, required=True,
        help="The desired full path to which the synthetic study object will be saved.."
    )
    parser.add_argument(
        '--read-depth', '-r', dest='read_depth',
        type=int, required=True,
        help="The overall read depth to simulate per timepoint."
    )
    parser.add_argument(
        '--a0', '-a0', dest='negbin_a0',
        type=float, required=True,
        help="The a0 parameter of the negbin parametrization."
    )
    parser.add_argument(
        '--a1', '-a1', dest='negbin_a1',
        type=float, required=True,
        help="The a1 parameter of the negbin parametrization."
    )
    parser.add_argument(
        '--qpcr-noise-scale', '-q', dest='qpcr_noise_scale',
        type=float, required=True,
        help="The qPCR noise scale (geometric noise stdev scaling)"
    )
    parser.add_argument(
        '--seed', dest='seed',
        type=int, required=True,
        help="The random seed to use for sampling randomness."
    )
    parser.add_argument('--sim-max', type=float, required=False, default=1e20)
    return parser.parse_args()


def main(
        study_path: Path,
        growth_path: Path,
        interactions_path: Path,
        perturbations_path: Path,
        coclusterings_path: Path,
        ground_truth_dir: Path,
        out_path: Path,
        replicate_out_path: Path,
        read_depth: int,
        negbin_a0: float,
        negbin_a1: float,
        qpcr_noise_scale: float,
        seed: int,
        sim_max: float,
):
    """
    Read the fixed-module inference, use the stored "filtered-state" (latent traj X) samples.

    :param fixed_module_pkl_path:
    :param read_depth:
    :return:
    """
    source_study = md2.Study.load(str(study_path))
    ground_truth_dir.mkdir(parents=True, exist_ok=True)
    print(f"Using ground truth dir {ground_truth_dir}, and then doing Discretized Euler integration.")

    glv_sim_path = ground_truth_dir / 'glv_best_sim.pkl'
    if glv_sim_path.exists():
        with open(glv_sim_path, 'rb') as f:
            glv_params, forward_sims = pkl.load(f)
    else:
        # Take the median parameter set from posterior.
        growths = np.load(growth_path)
        interactions = np.load(interactions_path)
        perturbations_map = np.load(perturbations_path)
        perturbations = [
            perturbations_map[pert.name]
            for pert in source_study.perturbations
        ]
        coclusterings = np.load(coclusterings_path)
        assert coclusterings.shape[0] == interactions.shape[0]  # ensure that the number of samples match

        glv_params, instance_coclusters, forward_sims = extract_glv_model(source_study, growths, interactions, perturbations, coclusterings, sim_max)
        plot_dir = ground_truth_dir / 'fwsim-plots'
        plot_dir.mkdir(exist_ok=True, parents=True)
        render_plots(source_study, forward_sims, plot_dir)
        with open(glv_sim_path, 'wb') as f:
            pkl.dump((glv_params, forward_sims), f)
        np.save(ground_truth_dir / "coclusters.npy", instance_coclusters)

        # Save the full forward simulation.
        for subj in source_study:
            sims, times = forward_simulate_subject(
                glv_params,
                source_study,
                subj,
                sim_max
            )
            np.savez(ground_truth_dir / f'forward_sim_full_{subj.name}.npz', sims=sims, times=times)

    rng = np.random.default_rng(seed)
    out_path.parent.mkdir(exist_ok=True, parents=True)
    synth_study = create_synthetic_dataset(
        source_study=source_study,
        forward_sims=forward_sims,
        sim_read_depth=read_depth,
        negbin_a0=negbin_a0,
        negbin_a1=negbin_a1,
        qpcr_noise_scale=qpcr_noise_scale,
        rng=rng,
    )

    synth_study.save(str(out_path))
    print(f"Saved semisynthetic dataset to {out_path}.")

    synth_repl_study = create_synthetic_replicates(
        source_study=source_study,
        forward_sims=forward_sims,
        source_subj_name='2',
        source_subj_timepoints=[8.0, 9.0, 10.0],  # timepoints
        num_physical_replicates=6,
        replicate_study_name='synthetic-replicate',
        negbin_a0=negbin_a0,
        negbin_a1=negbin_a1,
        num_reads=read_depth,
        qpcr_noise_scale=qpcr_noise_scale,
        rng=rng
    )
    replicate_out_path.parent.mkdir(exist_ok=True, parents=True)
    synth_repl_study.save(str(replicate_out_path))
    print(f"Saved replicates to {replicate_out_path}.")


if __name__ == "__main__":
    args = parse_args()
    main(
        study_path=Path(args.study_path),
        growth_path=Path(args.growth_path),
        interactions_path=Path(args.interactions_path),
        perturbations_path=Path(args.perts_path),
        coclusterings_path=Path(args.coclust_path),
        ground_truth_dir=Path(args.ground_truth_dir),
        out_path=Path(args.out_path),
        replicate_out_path=Path(args.replicate_out_path),
        read_depth=args.read_depth,
        negbin_a0=args.negbin_a0,
        negbin_a1=args.negbin_a1,
        qpcr_noise_scale=args.qpcr_noise_scale,
        seed=args.seed,
        sim_max=args.sim_max,
    )
