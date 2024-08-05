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
    parser.add_argument('-n', '--num_subjects', type=int, required=True,
                        help='<Required> The number of subjecs to simulate to lump into a single cohort.')

    parser.add_argument('--subsample_fwsim', type=int, required=False, default=0.1,
                        help='<Optional> Subsamples this fraction of poster samples from MDSINE(1 or 2).'
                             'Used for forward simulation-based metrics only.')
    return parser.parse_args()


# =============== Iterators through directory structures ==========
def result_dirs(results_base_dir: Path) -> Iterator[Tuple[int, int, str, Path]]:
    for read_depth, read_depth_dir in read_depth_dirs(results_base_dir):
        if read_depth == 1000:
            print("Skipping read depth 1000.")
            continue
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
    mat_file = result_dir / "mdsine1" / "BVS.mat"
    if not mat_file.exists():
        raise FileNotFoundError(str(mat_file))

    with h5py.File(mat_file, 'r') as f:
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
        return np.transpose(interactions, [0, 2, 1]), growths, indicator_probs.T


def mdsine2_output(result_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parser which extracts glv params from the specified results directory.
    :return:
    """
    mcmc = md2.BaseMCMC.load(str(result_dir / "mcmc.pkl"))
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

    rd = Ridge(glv.X, glv.T, glv.r_A, glv.r_B, glv.r_g, U, glv.scale)

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
            try:
                pred_growth, _ = regression_output(result_dir, _method, _regression_type)
                _add_entry(f'{_method}-{_regression_type}', _error_metric(pred_growth, true_growth))
            except FileNotFoundError:
                pass

        # MDSINE2 inference error eval
        try:
            _, growths, _ = mdsine2_output(result_dir / "mdsine2" / f"simulated-{noise_level}")
            errors = np.array([_error_metric(pred_growth, true_growth) for pred_growth in growths])
            _add_entry('MDSINE2', float(np.median(errors)))
        except FileNotFoundError:
            print("Skipping MDSINE2 (read depth={}, trial={}, noise={}), dir = {}".format(read_depth, trial_num, noise_level, result_dir))
            pass

        # MDSINE1 error
        try:
            _, growths, _ = mdsine1_output(result_dir)
            errors = np.array([_error_metric(pred_growth, true_growth) for pred_growth in growths])
            _add_entry('MDSINE1', float(np.median(errors)))
        except FileNotFoundError:
            pass

        # CLV inference error eval
        # _add_regression_entry("lra", "elastic_net")
        _add_regression_entry("glv", "elastic_net")
        _add_regression_entry("glv", "ridge")
        # _add_regression_entry("glv-ra", "elastic_net")
        # _add_regression_entry("glv-ra", "ridge")

    df = pd.DataFrame(df_entries)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    return df


def evaluate_interaction_strength_errors(true_interactions: np.ndarray, results_base_dir: Path) -> pd.DataFrame:
    df_entries = []

    def _error_metric(pred, truth) -> float:
        assert pred.shape[0] == pred.shape[1]  # Square matrix
        # return np.sqrt(np.mean(np.square(pred - truth)))

        # ===== Don't count diagonals (self interactions)
        np.fill_diagonal(pred, 0)
        np.fill_diagonal(truth, 0)
        num_entries = pred.shape[0] * (pred.shape[0] - 1)
        return np.sqrt((1 / num_entries) * np.sum(np.square(pred - truth)))

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
            try:
                _, pred_interaction = regression_output(result_dir, _method, _regression_type)
                _add_entry(f'{_method}-{_regression_type}', _error_metric(np.transpose(pred_interaction), true_interactions))
            except FileNotFoundError:
                pass

        # MDSINE2 inference error eval
        try:
            interactions, _, _ = mdsine2_output(result_dir / "mdsine2" / f"simulated-{noise_level}")
            errors = np.array([_error_metric(pred_interaction, true_interactions) for pred_interaction in interactions])
            _add_entry('MDSINE2', float(np.median(errors)))
            # pred_interaction = np.median(interactions, axis=0)
            # _add_entry('MDSINE2', _error_metric(pred_interaction, true_interactions))
        except FileNotFoundError:
            pass

        # MDSINE1 error
        try:
            interactions, _, _ = mdsine1_output(result_dir)
            errors = np.array([_error_metric(pred_interaction, true_interactions) for pred_interaction in interactions])
            _add_entry('MDSINE1', float(np.median(errors)))
        except FileNotFoundError:
            pass

        # CLV inference error eval
        # _add_regression_entry("lra", "elastic_net")
        _add_regression_entry("glv", "elastic_net")
        _add_regression_entry("glv", "ridge")
        # _add_regression_entry("glv-ra", "elastic_net")
        # _add_regression_entry("glv-ra", "ridge")

    df = pd.DataFrame(df_entries)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    return df


def evaluate_topology_errors(true_indicators: np.ndarray, results_base_dir: Path) -> pd.DataFrame:
    df_entries = []
    from sklearn.metrics import roc_auc_score

    def flatten_off_diagonal(A: np.ndarray) -> np.ndarray:
        return np.concatenate([
            A[np.triu_indices_from(A, k=1)],
            A[np.tril_indices_from(A, k=-1)]
        ])

    def _compute_auroc(_method: str, _p: np.ndarray):
        auc = roc_auc_score(
            flatten_off_diagonal(true_indicators),
            flatten_off_diagonal(_p),
        )
        return auc

    for read_depth, trial_num, noise_level, result_dir in result_dirs(results_base_dir):
        def _compute_roc_curve(_method: str, _pred: np.ndarray):
            """
            Computes a range of FPR/TPR pairs and adds them all to the dataframe.
            """
            df_entries.append({
                'Method': _method,
                'ReadDepth': read_depth,
                'Trial': trial_num,
                'NoiseLevel': noise_level,
                'AUROC': _compute_auroc(_method, _pred)
            })

        def _add_regression_entry(_method: str, _regression_type: str):
            try:
                interaction_p_values = regression_interaction_pvals(result_dir, _method, _regression_type)  # (N x N)
                # smalller p-value = more significant, so convert it into a "score" by passing 1-p.
                _compute_roc_curve(f'{_method}-{_regression_type}', 1.0 - interaction_p_values)
            except FileNotFoundError:
                pass

        # MDSINE2 inference error eval
        try:
            _, _, interaction_indicators = mdsine2_output(result_dir / 'mdsine2' / f'simulated-{noise_level}')
            indicator_probs = np.mean(interaction_indicators, axis=0)
            _compute_roc_curve('MDSINE2', indicator_probs)
        except FileNotFoundError:
            pass

        # MDSINE1 inference error eval
        try:
            _, _, indicator_probs = mdsine1_output(result_dir)
            _compute_roc_curve('MDSINE1', indicator_probs)
        except FileNotFoundError:
            pass

        # CLV inference error eval
        # Note: No obvious t-test implementation for elastic net regression.
        # _add_regression_entry("lra", "elastic_net")
        _add_regression_entry("glv", "elastic_net")
        _add_regression_entry("glv", "ridge")
        # _add_regression_entry("glv-ra", "elastic_net")
        # _add_regression_entry("glv-ra", "ridge")

    df = pd.DataFrame(df_entries)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    return df


def evaluate_fwsim_errors(true_growth: np.ndarray,
                          true_interactions: np.ndarray,
                          initial_min_value: float,
                          results_base_dir: Path,
                          dataset_dir: Path,
                          n_subjects: int,
                          subsample_frac: float) -> pd.DataFrame:
    """
    Generate a new subject and use it as a "holdout" dataset.
    :param true_growth:
    :param true_interactions:
    :param results_base_dir:
    :return:
    """
    def _error_metric(_pred_trajs, _true_trajs) -> float:
        _pred_trajs[_pred_trajs < 1.0] = 1.0
        _pred_trajs = np.log10(_pred_trajs)
        _true_trajs = np.log10(_true_trajs)
        return np.sqrt(np.mean(np.square(_pred_trajs - _true_trajs)))

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

        # load initial conditions.
        true_trajs = []
        initial_conds = []
        for subj_idx in range(n_subjects):
            subj_name = f'subj_{subj_idx}'
            initial_cond = np.load(dataset_dir / f'trial_{trial_num}' / f'{subj_name}.npz')['sims'][:, 0]
            initial_cond[initial_cond < initial_min_value] = initial_min_value
            initial_conds.append(initial_cond)

            true_traj_subj, _ = forward_sim(true_growth, true_interactions, initial_cond, dt=sim_dt, sim_max=sim_max, sim_t=sim_t)
            true_traj_subj = true_traj_subj[:, target_t_idx]
            true_trajs.append(true_traj_subj)
        true_trajs = np.stack(true_trajs, axis=0)  # (n_subj x n_taxa x n_times)

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
            pred_trajs = np.stack([
                np.median(
                    posterior_forward_sims(_pred_growths[subsample_idxs, :],
                                           _pred_interactions[subsample_idxs, :, :],
                                           initial_cond, sim_dt, sim_max, sim_t, t, target_t_idx),
                    axis=0
                )  # median trajectory
                for initial_cond in initial_conds
            ])
            _add_entry(_method, _error_metric(pred_trajs, true_trajs))

        def _eval_regression(_model_name: str, _regression_type: str):
            pkl_dir = result_dir / _model_name / _regression_type
            result_paths = list(pkl_dir.glob('*.pkl'))
            if len(result_paths) == 0:
                print(f"Unable to locate any .pkl files in {result_dir}.")
                return
            pkl_path = result_paths[0]
            pred_trajs = np.stack([
                np.transpose(regression_forward_simulate(pkl_path, init_abundance=initial_cond, times=target_t))
                for initial_cond in initial_conds
            ])
            _add_entry(f'{_model_name}-{_regression_type}', _error_metric(pred_trajs, true_trajs))

        # MDSINE2 error
        try:
            pred_interactions, pred_growths, _ = mdsine2_output(result_dir)
            _eval_mdsine('MDSINE2', pred_interactions, pred_growths)
        except FileNotFoundError:
            pass

        # MDSINE1 error
        try:
            pred_interactions, pred_growths, _ = mdsine1_output(result_dir)
            _eval_mdsine('MDSINE1', pred_interactions, pred_growths)
        except FileNotFoundError:
            pass

        # _eval_regression("lra", "elastic_net")
        _eval_regression("glv", "elastic_net")
        _eval_regression("glv", "ridge")
        # _eval_regression("glv-ra", "elastic_net")
        # _eval_regression("glv-ra", "ridge")
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
       forward simulates the gLV model using parameters inferred by the
       regression model

       init_abundance : N dimensional array containing the initial abundances
       times : T dimensional array containing the observation times

       @return
       T x N array containing the predicted abundances
    """
    def grad_fn(A, g, B, u):
        def fn(t, x):
            #for gLV we use log concetrations
            if B is None or u is None:
                return g + A.dot(np.exp(x))
            elif B is not None and u is not None:
                return g + A.dot(np.exp(x)) + B.dot(u)

        return fn

    with open(glv_pkl_loc, "rb") as f:
        glv = pickle.load(f)
    x_pred = np.zeros((times.shape[0], init_abundance.shape[0]))
    x_pred[0] = init_abundance
    xt = np.log(init_abundance)

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
        x_pred[t] = np.exp(xt)

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
    dyn = md2.model.gLVDynamicsSingleClustering(
        growth=None, interactions=None,
        perturbation_ends=[], perturbation_starts=[],
        start_day=0, sim_max=sim_max)
    initial_conditions = initial_conditions.reshape(-1, 1)

    dyn.growth = growth
    dyn.interactions = interactions
    dyn.perturbations = []

    x = md2.integrate(
        dynamics=dyn,
        initial_conditions=initial_conditions,
        dt=dt,
        final_day=sim_t,
        subsample=False
    )
    return x['X'], x['times']


def main():
    args = parse_args()
    growth, interactions, indicators = parse_ground_truth_params(Path(args.ground_truth_params))

    results_base_dir = Path(args.results_base_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    print(f"Outputs will be saved to {output_dir}.")

    print("Evaluating growth rate errors.")
    growth_rate_errors = evaluate_growth_rate_errors(growth, results_base_dir)
    out_path = output_dir / "growth_rate_errors.csv"
    growth_rate_errors.to_csv(output_dir / "growth_rate_errors.csv")
    print(f"Wrote growth rate errors to {out_path.name}.")

    print("Evaluating interaction strength errors.")
    interaction_strength_errors = evaluate_interaction_strength_errors(interactions, results_base_dir)
    out_path = output_dir / "interaction_strength_errors.csv"
    interaction_strength_errors.to_csv(output_dir / "interaction_strength_errors.csv")
    print(f"Wrote interaction strength errors to {out_path.name}.")

    print("Evaluating interaction topology errors.")
    topology_errors = evaluate_topology_errors(indicators, results_base_dir)
    out_path = output_dir / "topology_errors.csv"
    topology_errors.to_csv(output_dir / "topology_errors.csv")
    print(f"Wrote interaction topology errors to {out_path.name}.")

    print("Evaluating forward-simulation errors.")
    dataset_dir = Path("/data/local/youn/MDSINE2_Paper/datasets/synthetic/data")
    holdout_trajectory_errors = evaluate_fwsim_errors(
        growth,
        interactions,
        100.0,
        results_base_dir,
        dataset_dir,
        args.num_subjects,
        args.subsample_fwsim
    )
    out_path = output_dir / "holdout_trajectory_errors.csv"
    holdout_trajectory_errors.to_csv(out_path)
    print(f"Wrote heldout trajectory prediction errors to {out_path.name}.")


if __name__ == "__main__":
    main()
