from typing import Tuple, Iterator
from pathlib import Path
import argparse
import pickle

import numpy as np
import pandas as pd
import mdsine2 as md2
from mdsine2.names import STRNAMES


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--ground_truth_params', type=str, required=True,
                        help='<Required> The path to a .npz file containing `growth_rates`, `interactions` arrays')
    return parser.parse_args()


# =============== Iterators through directory structures ==========
def result_dirs(results_base_dir: Path) -> Iterator[Tuple[int, int, str, Path]]:
    for read_depth, read_depth_dir in read_depth_dirs(results_base_dir):
        for trial_num, trial_dir in trial_dirs(read_depth_dir):
            for noise_level, noise_level_dir in noise_level_dirs(trial_dir):
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
def mdsine_output(result_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parser which extracts glv params from the specified results directory.
    :return:
    """
    mcmc = md2.BaseMCMC.load(str(result_dir / "mdsine2" / "mcmc.pkl"))
    interactions = mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk(section='posterior')
    growths = mcmc.graph[STRNAMES.GROWTH_VALUE].get_trace_from_disk(section='posterior')
    interaction_indicators = mcmc.graph[STRNAMES.CLUSTER_INTERACTION_INDICATOR].get_trace_from_disk(section='posterior')
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


def regression_interaction_pvals(result_dir: Path, model_name: str, regression_type: str) -> np.ndarray:
    """
    Generator which yields growth/interaction params estimated through regression.
    :param result_dir:
    :param model_name:
    :param regression_type:
    :return:
    """
    raise NotImplementedError()


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
        _, growths, _ = mdsine_output(result_dir)
        pred_growth = np.median(growths, axis=0)
        _add_entry('MDSINE2', _error_metric(pred_growth, true_growth))

        # CLV inference error eval
        _add_regression_entry("lra", "elastic_net")
        _add_regression_entry("glv", "elastic_net")
        _add_regression_entry("glv", "ridge")
        _add_regression_entry("glv-ra", "elastic_net")
        _add_regression_entry("glv-ra", "ridge")

    return pd.DataFrame(df_entries)


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
        interactions, _, _ = mdsine_output(result_dir)
        pred_interaction = np.median(interactions, axis=0)
        _add_entry('MDSINE2', _error_metric(pred_interaction, true_interactions))

        # CLV inference error eval
        _add_regression_entry("lra", "elastic_net")
        _add_regression_entry("glv", "elastic_net")
        _add_regression_entry("glv", "ridge")
        _add_regression_entry("glv-ra", "elastic_net")
        _add_regression_entry("glv-ra", "ridge")

    return pd.DataFrame(df_entries)


def evaluate_topology_errors(true_indicators: np.ndarray, results_base_dir: Path) -> pd.DataFrame:
    df_entries = []

    def _false_positive_rate(pred, truth) -> float:
        return np.sum(pred & np.logical_not(truth)) / len(truth)

    def _true_positive_rate(pred, truth) -> float:
        raise np.sum(pred & truth) / len(truth)

    for read_depth, trial_num, noise_level, result_dir in result_dirs(results_base_dir):
        def _compute_roc_curve(_method: str, _pvalues: np.ndarray):
            """
            Computes a range of FPR/TPR pairs and adds them all to the dataframe.
            :param _pvalues: An (N_taxa x N_taxa) matrix of p-values for each pair of interactions.
            """
            q = np.linspace(0, 1, 1000)
            preds = np.expand_dims(_pvalues, axis=2) > np.expand_dims(q, axis=(0, 1))
            for i in range(len(q)):
                df_entries.append({
                    'Method': _method,
                    'ReadDepth': read_depth,
                    'Trial': trial_num,
                    'NoiseLevel': noise_level,
                    'FPR': _false_positive_rate(preds[:, :, i], true_indicators),
                    'TPR': _true_positive_rate(preds[:, :, i], true_indicators)
                })

        def _add_regression_entry(_method: str, _regression_type: str):
            interaction_p_values = regression_interaction_pvals(result_dir, _method, _regression_type)  # (N x N)
            _compute_roc_curve(f'{_method}-{_regression_type}', interaction_p_values)

        # MDSINE2 inference error eval
        _, _, interaction_indicators = mdsine_output(result_dir)
        indicator_pvals = np.mean(interaction_indicators, axis=0)
        _compute_roc_curve('MDSINE2', indicator_pvals)

        # CLV inference error eval
        _add_regression_entry("lra", "elastic_net")
        _add_regression_entry("glv", "elastic_net")
        _add_regression_entry("glv", "ridge")
        _add_regression_entry("glv-ra", "elastic_net")
        _add_regression_entry("glv-ra", "ridge")

    return pd.DataFrame(df_entries)


def evaluate_holdout_trajectory_errors(true_growth: np.ndarray, true_interactions: np.ndarray, results_base_dir: Path) -> pd.DataFrame:
    """
    Generate a new subject and use it as a "holdout" dataset.
    :param true_growth:
    :param true_interactions:
    :param results_base_dir:
    :return:
    """
    raise NotImplementedError()


def parse_taxa_list(taxa_path: Path):
    taxon = []
    with open(taxa_path, 'r') as f:
        for line in f:
            line = line.rstrip()
            if len(line) > 0:
                taxon.append(line)
    return taxon


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

    holdout_trajectory_errors = evaluate_holdout_trajectory_errors(growth, interactions, results_base_dir)
    holdout_trajectory_errors.to_csv(output_dir / "holdout_trajectory_errors.csv")
    print(f"Wrote heldout trajectory prediction errors.")


if __name__ == "__main__":
    main()
