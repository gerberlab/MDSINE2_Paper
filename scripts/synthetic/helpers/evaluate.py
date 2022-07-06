from typing import Tuple, Iterator
from pathlib import Path
import argparse
import pickle

import numpy as np
import h5py
import pandas as pd
import scipy.stats
from scipy.integrate import solve_ivp

import mdsine2 as md2
from mdsine2.names import STRNAMES
from regression_analyzer import Ridge
from generalized_lotka_volterra import GeneralizedLotkaVolterra


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--results_base_dir', type=str, required=True,
                        help='<Required> The path to the output/ base directory containing the results of all '
                             'inference runs.')
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help='<Required> The desired path to which all evaluation metrics will be saved.')
    parser.add_argument('-g', '--ground_truth_params', type=str, required=True,
                        help='<Required> The path to a .npz file containing `growth_rates`, `interactions` arrays')
    parser.add_argument('-m', '--initial_cond_mean', type=float, required=True,
                        help='<Required> The mean of the initial condition distribution.')
    parser.add_argument('-s', '--initial_cond_std', type=float, required=True,
                        help='<Required> The standard deviation of the initial condition distribution.')

    parser.add_argument('--subsample_fwsim', type=int, required=False, default=0.1,
                        help='<Optional> Subsamples this fraction of poster samples from MDSINE(1 or 2).'
                             'Used for forward simulation-based metrics only.')
    return parser.parse_args()


# =============== Iterators through directory structures ==========
def result_dirs(results_base_dir: Path) -> Iterator[Tuple[int, int, str, Path]]:
    for read_depth, read_depth_dir in read_depth_dirs(results_base_dir):
        for trial_num, trial_dir in trial_dirs(read_depth_dir):
            for noise_level, noise_level_dir in noise_level_dirs(trial_dir):
                print(f"Yielding (read_depth: {read_depth}, trial: {trial_num}, noise: {noise_level})")
                yield read_depth, trial_num, noise_level, noise_level_dir


def read_depth_dirs(results_base_dir: Path) -> Iterator[Tuple[int, Path]]:
    for child_dir in results_base_dir.glob("reads_*"):
        if not child_dir.is_dir():
            raise RuntimeError(f"Expected child `{child_dir}` to be a directory.")

        read_depth = int(child_dir.name.split("_")[1])
        yield read_depth, child_dir


def trial_dirs(read_depth_dir: Path) -> Iterator[Tuple[int, Path]]:
    for child_dir in read_depth_dir.glob("trial_*"):
        if not child_dir.is_dir():
            raise RuntimeError(f"Expected child `{child_dir}` to be a directory.")

        trial_num = int(child_dir.name.split("_")[1])
        yield trial_num, child_dir


def noise_level_dirs(trial_dir: Path) -> Iterator[Tuple[str, Path]]:
    for child_dir in trial_dir.glob("*_noise"):
        if not child_dir.is_dir():
            raise RuntimeError(f"Expected child `{child_dir}` to be a directory.")

        noise_level = child_dir.name.split("_")[0]
        yield noise_level, child_dir


# ========================= Method output iterators =======================
def mdsine1_output(result_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    :param result_dir:
    :return: posterior mean interactions, posterior mean growth, posterior mean indicator.
    """
    with h5py.File(result_dir / "mdsine1" / "BVS.mat", 'r') as f:
        theta_mean = np.array(f['Theta_select'])
        theta_samples = f['Theta_samples_select'][0]

        all_species = [
            ''.join(chr(x) for x in f[ref])
            for ref in f['species_names'][0]
        ]
        filtered_species = [
            ''.join(chr(x) for x in f[ref])
            for ref in f['species_names_filtered_total'][0]
        ]
        filtered_indices = [all_species.index(sp) for sp in filtered_species]

        n_samples = len(theta_samples)
        n_taxa = len(all_species)

        growths = np.zeros(shape=(n_samples, n_taxa), dtype=float)
        interactions = np.zeros(shape=(n_samples, n_taxa, n_taxa), dtype=float)
        indicator_probs = np.zeros(shape=(n_taxa, n_taxa))
        indicator_probs[np.ix_(filtered_indices, filtered_indices)] = np.array(f['Theta_select_probs'])

        for n in range(n_samples):
            ref_n = theta_samples[n]
            theta_n = np.array(f[ref_n])

            growths[n, filtered_indices] = theta_n[0, :]
            interaction_slice = interactions[n]
            interaction_slice[np.ix_(filtered_indices, filtered_indices)] = theta_n[1:, :]
        return interactions, growths, indicator_probs


def mdsine2_output(result_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parser which extracts glv params from the specified results directory.
    :return:
    """
    mcmc = md2.BaseMCMC.load(str(result_dir / "mdsine2" / "mcmc.pkl"))
    growths = mcmc.graph[STRNAMES.GROWTH_VALUE].get_trace_from_disk(section='posterior')
    interaction_indicators = mcmc.graph[STRNAMES.CLUSTER_INTERACTION_INDICATOR].get_trace_from_disk(section='posterior')
    self_interactions = mcmc.graph[STRNAMES.SELF_INTERACTION_VALUE].get_trace_from_disk(section='posterior')
    interactions = mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk(section='posterior')

    self_interactions = -np.absolute(self_interactions)
    interactions[np.isnan(interactions)] = 0
    for interaction, indicator in zip(interactions, interaction_indicators):
        interaction[indicator == 0] = 0.
    for i in range(self_interactions.shape[1]):
        interactions[:, i, i] = self_interactions[:, i]
    return interactions, growths, interaction_indicators


def regression_output(result_dir: Path, model_name: str, regression_type: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generator which yields growth/interaction params estimated through regression.
    :param result_dir: The directory containing the particular regressions output pickle file(s).
    :param model_name:
    :param regression_type:
    :return:
    """
    result_dir = result_dir / model_name / regression_type
    result_paths = list(result_dir.glob('*.pkl'))
    if len(result_paths) == 0:
        raise FileNotFoundError(f"Unable to locate any .pkl files in {result_dir}.")
    result_path = result_paths[0]

    with open(result_path, 'rb') as f:
        res = pickle.load(f)
    est_interactions = res.A
    est_growth = res.g
    return est_growth, est_interactions


# noinspection PyPep8Naming
def regression_interaction_pvals(result_dir: Path, model_name: str, regression_type: str) -> np.ndarray:
    """
    Generator which yields growth/interaction params estimated through regression.
    :param result_dir:
    :param model_name:
    :param regression_type:
    :return:
    """
    result_dir = result_dir / model_name / regression_type
    result_paths = list(result_dir.glob('*.pkl'))
    if len(result_paths) == 0:
        raise FileNotFoundError(f"Unable to locate any .pkl files in {result_dir}.")
    result_path = result_paths[0]

    with open(result_path, "rb") as f:
        glv: GeneralizedLotkaVolterra = pickle.load(f)

    U = glv.U
    if np.sum(np.vstack(U)) == 0:
        U = None

    rd = Ridge(glv.X, glv.T, glv.r_A, glv.r_B, glv.r_g, U)

    """ p_interaction[i, j]: p-value associated with the effect of bug i on bug j """
    p_interaction, p_growth, p_perturbation = rd.significance_test()
    return p_interaction


# ======================== Error metrics =======================
def parse_ground_truth_params(params_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    params = np.load(str(params_path))
    growth = params['growth_rates']
    interactions = params['interactions']
    indicators = (interactions != 0.0)
    return growth, interactions, indicators


def evaluate_growth_rate_errors(true_growth: np.ndarray, results_base_dir: Path) -> pd.DataFrame:
    df_entries = []

    def _error_metric(pred, truth) -> float:
        return np.sqrt(np.mean(np.square(pred - truth)))

    for read_depth, trial_num, noise_level, result_dir in result_dirs(results_base_dir):
        def _add_entry(_method: str, _err: float):
            df_entries.append({
                'Method': _method,
                'ReadDepth': read_depth,
                'Trial': trial_num,
                'NoiseLevel': noise_level,
                'Error': _err
            })

        def _add_regression_entry(_method: str, _regression_type: str):
            pred_growth, _ = regression_output(result_dir, _method, _regression_type)
            _add_entry(f'{_method}-{_regression_type}', _error_metric(pred_growth, true_growth))

        # MDSINE2 inference error eval
        _, growths, _ = mdsine2_output(result_dir)
        errors = np.array([_error_metric(pred_growth, true_growth) for pred_growth in growths])
        _add_entry('MDSINE2', float(np.median(errors)))

        # MDSINE1 error
        _, growths, _ = mdsine1_output(result_dir)
        errors = np.array([_error_metric(pred_growth, true_growth) for pred_growth in growths])
        _add_entry('MDSINE1', float(np.median(errors)))

        # CLV inference error eval
        _add_regression_entry("lra", "elastic_net")
        _add_regression_entry("glv", "elastic_net")
        _add_regression_entry("glv", "ridge")
        _add_regression_entry("glv-ra", "elastic_net")
        _add_regression_entry("glv-ra", "ridge")

    df = pd.DataFrame(df_entries)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    return df


def evaluate_interaction_strength_errors(true_interactions: np.ndarray, results_base_dir: Path) -> pd.DataFrame:
    df_entries = []

    def _error_metric(pred, truth) -> float:
        return np.sqrt(np.mean(np.square(pred - truth)))

    for read_depth, trial_num, noise_level, result_dir in result_dirs(results_base_dir):
        def _add_entry(_method: str, _err: float):
            df_entries.append({
                'Method': _method,
                'ReadDepth': read_depth,
                'Trial': trial_num,
                'NoiseLevel': noise_level,
                'Error': _err
            })

        def _add_regression_entry(_method: str, _regression_type: str):
            _, pred_interaction = regression_output(result_dir, _method, _regression_type)
            _add_entry(f'{_method}-{_regression_type}', _error_metric(np.transpose(pred_interaction), true_interactions))

        # MDSINE2 inference error eval
        interactions, _, _ = mdsine2_output(result_dir)
        errors = np.array([_error_metric(pred_interaction, true_interactions) for pred_interaction in interactions])
        _add_entry('MDSINE2', float(np.median(errors)))
        # pred_interaction = np.median(interactions, axis=0)
        # _add_entry('MDSINE2', _error_metric(pred_interaction, true_interactions))

        # MDSINE1 error
        interactions, _, _ = mdsine1_output(result_dir)
        errors = np.array([_error_metric(pred_interaction, true_interactions) for pred_interaction in interactions])
        _add_entry('MDSINE1', float(np.median(errors)))

        # CLV inference error eval
        _add_regression_entry("lra", "elastic_net")
        _add_regression_entry("glv", "elastic_net")
        _add_regression_entry("glv", "ridge")
        _add_regression_entry("glv-ra", "elastic_net")
        _add_regression_entry("glv-ra", "ridge")

    df = pd.DataFrame(df_entries)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    return df


def evaluate_topology_errors(true_indicators: np.ndarray, results_base_dir: Path) -> pd.DataFrame:
    df_entries = []

    def _false_positive_rate(pred, truth) -> float:
        return 1 - _true_negative_rate(pred, truth)

    def _true_negative_rate(pred, truth) -> float:
        return np.sum(~pred & ~truth) / np.sum(~truth)

    def _true_positive_rate(pred, truth) -> float:
        return np.sum(pred & truth) / np.sum(truth)

    for read_depth, trial_num, noise_level, result_dir in result_dirs(results_base_dir):
        def _compute_roc_curve(_method: str, _p: np.ndarray, _q: np.ndarray, use_greater_than: bool):
            """
            Computes a range of FPR/TPR pairs and adds them all to the dataframe.
            :param _p: An (N_taxa x N_taxa) matrix of p-values for each pair of interactions.
            """
            if use_greater_than:
                preds = np.expand_dims(_p, axis=2) > np.expand_dims(_q, axis=(0, 1))
            else:
                preds = np.expand_dims(_p, axis=2) < np.expand_dims(_q, axis=(0, 1))
            for i in range(len(_q)):
                preds_i = preds[:, :, i]
                np.fill_diagonal(true_indicators, 0)
                np.fill_diagonal(preds_i, 0)
                df_entries.append({
                    'Method': _method,
                    'ReadDepth': read_depth,
                    'Trial': trial_num,
                    'NoiseLevel': noise_level,
                    'FPR': _false_positive_rate(preds_i, true_indicators),
                    'TPR': _true_positive_rate(preds_i, true_indicators)
                })

        def _add_regression_entry(_method: str, _regression_type: str):
            interaction_p_values = regression_interaction_pvals(result_dir, _method, _regression_type)  # (N x N)
            _compute_roc_curve(f'{_method}-{_regression_type}', np.log(interaction_p_values), np.linspace(-100., 0., 1000), use_greater_than=False)

        # MDSINE2 inference error eval
        _, _, interaction_indicators = mdsine2_output(result_dir)
        indicator_pvals = np.mean(interaction_indicators, axis=0)
        _compute_roc_curve('MDSINE2', indicator_pvals, np.linspace(0., 1., 1000), use_greater_than=True)

        # MDSINE1 inference error eval
        _, _, indicator_probs = mdsine1_output(result_dir)
        _compute_roc_curve('MDSINE1', indicator_probs, np.linspace(0., 1., 1000), use_greater_than=True)

        # CLV inference error eval
        # Note: No obvious t-test implementation for elastic net regression.
        # _add_regression_entry("lra", "elastic_net")
        # _add_regression_entry("glv", "elastic_net")
        _add_regression_entry("glv", "ridge")
        # _add_regression_entry("glv-ra", "elastic_net")
        _add_regression_entry("glv-ra", "ridge")

    df = pd.DataFrame(df_entries)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    return df


def evaluate_holdout_trajectory_errors(true_growth: np.ndarray,
                                       true_interactions: np.ndarray,
                                       init_rv: scipy.stats.rv_continuous,
                                       results_base_dir: Path,
                                       subsample_frac: float) -> pd.DataFrame:
    """
    Generate a new subject and use it as a "holdout" dataset.
    :param true_growth:
    :param true_interactions:
    :param init_rv: The (1-d) distribution from which the initial condition will be sampled (iid for each taxon)
    :param results_base_dir:
    :return:
    """
    def _error_metric(_pred_traj, _true_traj) -> float:
        return np.sqrt(np.mean(np.square(_pred_traj - _true_traj)))

    """ Simulation parameters """
    sim_seed = 0
    sim_dt = 0.01
    sim_max = 1e20
    sim_t = 20
    t = np.arange(0., sim_t, sim_dt)
    target_t_idx = np.arange(len(t) // 2, len(t), int(1 / sim_dt))
    target_t = t[target_t_idx]

    """ Extraction of errors/simulations """
    df_entries = []
    for read_depth, trial_num, noise_level, result_dir in result_dirs(results_base_dir):
        sim_seed += 1
        np.random.seed(sim_seed)
        initial_cond = init_rv.rvs(size=len(true_growth))
        true_traj, _ = forward_sim(true_growth, true_interactions, initial_cond, dt=sim_dt, sim_max=sim_max, sim_t=sim_t)
        true_traj = true_traj[:, target_t_idx]

        def _add_entry(_method: str, _err: float):
            df_entries.append({
                'Method': _method,
                'ReadDepth': read_depth,
                'Trial': trial_num,
                'NoiseLevel': noise_level,
                'Error': _err
            })

        def _eval_mdsine(_method: str, _pred_interactions: np.ndarray, _pred_growths: np.ndarray):
            subsample_idxs = np.arange(0, _pred_interactions.shape[0], int(subsample_frac * _pred_interactions.shape[0]))
            pred_traj = np.median(
                posterior_forward_sims(_pred_growths[subsample_idxs, :],
                                       _pred_interactions[subsample_idxs, :, :],
                                       initial_cond, sim_dt, sim_max, sim_t, t, target_t_idx),
                axis=0
            )
            _add_entry(_method, _error_metric(pred_traj, true_traj))

        def _eval_regression(_model_name: str, _regression_type: str):
            pkl_dir = result_dir / _model_name / _regression_type
            result_paths = list(pkl_dir.glob('*.pkl'))
            if len(result_paths) == 0:
                raise FileNotFoundError(f"Unable to locate any .pkl files in {result_dir}.")
            pkl_path = result_paths[0]
            pred_traj = np.transpose(regression_forward_simulate(pkl_path, init_abundance=initial_cond, times=target_t))
            _add_entry(f'{_model_name}-{_regression_type}', _error_metric(pred_traj, true_traj))

        # MDSINE2 error
        pred_interactions, pred_growths, _ = mdsine2_output(result_dir)
        _eval_mdsine('MDSINE2', pred_interactions, pred_growths)

        # MDSINE1 error
        #pred_interactions, pred_growths, _ = mdsine1_output(result_dir)
        #_eval_mdsine('MDSINE1', pred_interactions, pred_growths)

        _eval_regression("lra", "elastic_net")
        _eval_regression("glv", "elastic_net")
        _eval_regression("glv", "ridge")
        _eval_regression("glv-ra", "elastic_net")
        _eval_regression("glv-ra", "ridge")
    df = pd.DataFrame(df_entries)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    return df


def regression_forward_simulate(glv_pkl_loc: Path,
                                init_abundance: np.ndarray,
                                times: np.ndarray,
                                perturbation_info: np.ndarray = None) -> np.ndarray:
    """
       forward simulates the trajectory using parameters inferred by the
       regression model

       init_abundance : N dimensional array containing the initial abundances
       times : T dimensional array containing the observation times

       @return
       T x N array containing the predicted abundances
    """
    def grad_fn(A, g, B, u):
        def fn(t, x):
            if B is None or u is None:
                return g + A.dot(x)
            elif B is not None and u is not None:
                return g + A.dot(x) + B.dot(u)

        return fn

    with open(glv_pkl_loc, "rb") as f:
        glv = pickle.load(f)
    x_pred = np.zeros((times.shape[0], init_abundance.shape[0]))
    x_pred[0] = init_abundance
    xt = init_abundance

    A, B, g = glv.A, glv.B, glv.g

    # no valid perturbation effects
    if np.sum(B) == 0:
        B = None

    for t in range(1, times.shape[0]):
        if perturbation_info is not None:
            grad = grad_fn(A, g, B, perturbation_info[t-1])
        else:
            grad = grad_fn(A, g, None, None)
        dt = times[t] - times[t-1]
        ivp = solve_ivp(grad, (0,0+dt), xt, method="RK45")
        xt = ivp.y[:,-1]
        x_pred[t] = xt

    return x_pred


def posterior_forward_sims(growths, interactions, initial_conditions, dt, sim_max, sim_t, expected_t, target_time_idxs) -> np.ndarray:
    """
    :param growths:
    :param interactions:
    :param initial_conditions:
    :param dt:
    :param sim_max:
    :param sim_t:
    :param expected_t: The array of timepoints we expect to find. Used for validation after forward simulation.
    :param target_time_idxs: An array of indexes of timepoints to be extracted.
    :return: M x N x T array of simulated trajs (M = # of MCMC samples, N = # of taxa, T = # timepoints)
    """
    fwsims = np.empty(
        shape=(growths.shape[0], growths.shape[-1], len(target_time_idxs)),
        dtype=float
    )
    for gibbs_idx in range(growths.shape[0]):
        growth = growths[gibbs_idx]
        interaction = interactions[gibbs_idx]
        _x, _t = forward_sim(growth, interaction, initial_conditions, dt, sim_max, sim_t)
        assert len(expected_t) <= len(_t)

        fwsims[gibbs_idx, :, :] = _x[:, target_time_idxs]

    # print(
    #     np.sum(np.isnan(fwsims), axis=0)
    # )
    return fwsims


def forward_sim(growth,
                interactions,
                initial_conditions,
                dt,
                sim_max,
                sim_t) -> Tuple[np.ndarray, np.ndarray]:
    """
    Forward simulate with the given dynamics. First start with the perturbation
    off, then on, then off.

    Parameters
    ----------
    growth : np.ndarray(n_gibbs, n_taxa)
        Growth parameters
    interactions : np.ndarray(n_gibbs, n_taxa, n_taxa)
        Interaction parameters
    initial_conditions : np.ndarray(n_taxa)
        Initial conditions of the taxa
    dt : float
        Step size to forward simulate with
    sim_max : float, None
        Maximum clip for forward sim
    sim_t : float
        Total number of days

    :return: N x T array of simulated trajs (N = # of taxa, T = # timepoints)
    """
    dyn = md2.model.gLVDynamicsSingleClustering(growth=None, interactions=None,
                                                perturbation_ends=[], perturbation_starts=[],
                                                start_day=0, sim_max=sim_max)
    initial_conditions = initial_conditions.reshape(-1, 1)

    dyn.growth = growth
    dyn.interactions = interactions
    dyn.perturbations = []

    x = md2.integrate(dynamics=dyn, initial_conditions=initial_conditions,
                      dt=dt, n_days=sim_t, subsample=False)
    return x['X'], x['times']


def main():
    args = parse_args()
    growth, interactions, indicators = parse_ground_truth_params(Path(args.ground_truth_params))

    results_base_dir = Path(args.results_base_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    print(f"Outputs will be saved to {output_dir}.")

    growth_rate_errors = evaluate_growth_rate_errors(growth, results_base_dir)
    growth_rate_errors.to_csv(output_dir / "growth_rate_errors.csv")
    print(f"Wrote growth rate errors.")

    interaction_strength_errors = evaluate_interaction_strength_errors(interactions, results_base_dir)
    interaction_strength_errors.to_csv(output_dir / "interaction_strength_errors.csv")
    print(f"Wrote interaction strength errors.")

    topology_errors = evaluate_topology_errors(indicators, results_base_dir)
    topology_errors.to_csv(output_dir / "topology_errors.csv")
    print(f"Wrote interaction topology errors.")

    init_dist = scipy.stats.norm(loc=args.initial_cond_mean, scale=args.initial_cond_std)
    holdout_trajectory_errors = evaluate_holdout_trajectory_errors(growth, interactions,
                                                                   init_dist, results_base_dir, args.subsample_fwsim)
    holdout_trajectory_errors.to_csv(output_dir / "holdout_trajectory_errors.csv")
    print(f"Wrote heldout trajectory prediction errors.")


if __name__ == "__main__":
    main()
