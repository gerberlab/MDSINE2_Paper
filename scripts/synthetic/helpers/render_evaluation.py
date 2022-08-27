from pathlib import Path
import argparse
from typing import List, Dict

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sb
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help='<Required> The target output directory. Should contain the results of `evaluate` script.')
    parser.add_argument('-f', '--format', type=str, required=False,
                        default='pdf',
                        help='<Optional> The plot image format (Default: pdf)')

    parser.add_argument('-fw', '--figure_width', type=int, required=False,
                        default=25,
                        help='<Optional> The figure width (Default: 16)')
    parser.add_argument('-fh', '--figure_height', type=int, required=False,
                        default=5,
                        help='<Optional> The figure height (Default: 5)')
    parser.add_argument('-r', '--read_depth', type=int, required=False,
                        default=25000,
                        help='<Optional> The read depth to render (Default: 25000)')
    return parser.parse_args()


def load_df(df_path: Path) -> pd.DataFrame:
    df = pd.read_csv(df_path)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    df = df.sort_values(['ReadDepth', 'NoiseLevel'], ascending=[True, True])
    return df


def render_growth_rate_errors(df: pd.DataFrame, ax, order, palette):
    df.loc[:, 'x'] = df['ReadDepth'].astype(str) + ' reads\n' + df['NoiseLevel'].astype(str) + ' noise'
    if df.shape[0] > 0:
        sb.barplot(
            data=df,
            x='x',
            y='Error',
            hue='Method',
            ax=ax,
            hue_order=order,
            palette=palette
        )
    ax.set_ylabel('RMSE')
    ax.set_xlabel(None)
    ax.set_yscale('log')
    ax.set_title('Growth Rates')
    ax.legend([], [], frameon=False)


def render_interaction_strength_errors(df: pd.DataFrame, ax, order, palette):
    df.loc[:, 'x'] = df['ReadDepth'].astype(str) + ' reads\n' + df['NoiseLevel'].astype(str) + ' noise'
    if df.shape[0] > 0:
        sb.barplot(
            data=df,
            x='x',
            y='Error',
            hue='Method',
            ax=ax,
            hue_order=order,
            palette=palette
        )
    ax.set_ylabel('RMSE')
    ax.set_xlabel(None)
    ax.set_yscale('log')
    ax.set_title('Interaction Strengths')
    ax.legend([], [], frameon=False)


def render_holdout_trajectory_errors(df: pd.DataFrame, ax, order, palette):
    df.loc[:, 'x'] = df['ReadDepth'].astype(str) + ' reads\n' + df['NoiseLevel'].astype(str) + ' noise'
    if df.shape[0] > 0:
        sb.barplot(
            data=df,
            x='x',
            y='Error',
            hue='Method',
            ax=ax,
            hue_order=order,
            palette=palette
        )
    ax.set_ylabel('RMSE')
    ax.set_xlabel(None)
    ax.set_title('Holdout Trajectories')
    ax.legend([], [], frameon=False)


def render_topology_errors(df: pd.DataFrame, ax, order, palette):
    def auroc(_df):
        _df = _df[['FPR', 'TPR']].groupby('FPR').max().reset_index()
        _df = _df.sort_values('FPR', ascending=True)
        fpr = np.concatenate([[0.], _df['FPR'].to_numpy(), [1.]])
        tpr = np.concatenate([[0.], _df['TPR'].to_numpy(), [1.]])
        return scipy.integrate.trapz(
            y=tpr,
            x=fpr,
        )

    area_df = df.groupby(['Method', 'ReadDepth', 'Trial', 'NoiseLevel']).apply(auroc).rename('Error').reset_index()
    area_df.loc[:, 'x'] = area_df['ReadDepth'].astype(str) + ' reads\n' + area_df['NoiseLevel'].astype(str) + ' noise'

    if area_df.shape[0] > 0:
        sb.barplot(
            data=area_df,
            x='x',
            y='Error',
            hue='Method',
            ax=ax,
            hue_order=order,
            palette=palette
        )
    ax.set_ylabel('AUC-ROC')
    ax.set_xlabel(None)
    ax.set_title('Network Structure')
    ax.legend([], [], frameon=False)


def render_all(fig: plt.Figure, dataframe_dir: Path, method_order: List[str], palette: Dict[str, np.ndarray], read_depth: int):

    # fig, axes = plt.subplots(1, 4, figsize=figsize)
    # fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    # spec = fig.add_gridspec(ncols=4, nrows=4, height_ratios=[1, 10, 1, 10], width_ratios=[1, 1, 1, 1])

    # ax0 = fig.add_subplot(spec[0, :2])
    # ax1, ax2 = fig.add_subplot(spec[1, 0]), fig.add_subplot(spec[1, 1])

    spec = fig.add_gridspec(ncols=5, nrows=1, width_ratios=[1, 1, 1, 1, 0.08], wspace=0.6)
    axes = [fig.add_subplot(spec[0, i]) for i in range(4)]
    df = load_df(dataframe_dir / "growth_rate_errors.csv")
    section = df.loc[df['ReadDepth'] == read_depth]
    render_growth_rate_errors(
        section,
        ax=axes[0],
        order=method_order,
        palette=palette
    )

    # ax0 = fig.add_subplot(spec[0, 2:])
    # ax1, ax2 = fig.add_subplot(spec[1, 2]), fig.add_subplot(spec[1, 3])
    df = load_df(dataframe_dir / "interaction_strength_errors.csv")
    section = df.loc[df['ReadDepth'] == read_depth]
    render_interaction_strength_errors(
        section,
        ax=axes[1],
        order=method_order,
        palette=palette
    )

    # ax0 = fig.add_subplot(spec[2, :2])
    # ax1, ax2 = fig.add_subplot(spec[3, 0]), fig.add_subplot(spec[3, 1])
    df = load_df(dataframe_dir / 'holdout_trajectory_errors.csv')
    section = df.loc[df['ReadDepth'] == read_depth]
    render_holdout_trajectory_errors(
        section,
        ax=axes[2],
        order=method_order,
        palette=palette
    )

    # ax0 = fig.add_subplot(spec[2, 2:])
    # ax1, ax2 = fig.add_subplot(spec[3, 2]), fig.add_subplot(spec[3, 3])

    topology_method_order = ['MDSINE2', 'MDSINE1', 'glv-ridge', 'glv-ra-ridge']
    df = load_df(dataframe_dir / "topology_errors.csv")
    section = df.loc[df['ReadDepth'] == read_depth]
    render_topology_errors(
        section,
        ax=axes[3],
        order=topology_method_order,
        palette=palette
    )


def draw_legend(fig, method_order: List[str], palette: Dict[str, np.ndarray]):
    legend_elements = [
        Patch(facecolor=palette[m], edgecolor=palette[m], label=m)
        for m in method_order
    ]

    fig.legend(
        handles=legend_elements,
        labels=method_order,
        bbox_to_anchor=(0.88, 0.5),
        markerfirst=False,
        loc='lower left'
    )


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)

    read_depth = args.read_depth

    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    method_order = ['MDSINE2', 'MDSINE1', 'lra-elastic_net', 'glv-elastic_net', 'glv-ridge', 'glv-ra-elastic_net', 'glv-ra-ridge']
    color_palette = {method: c for method, c in zip(method_order, default_colors)}

    fig = plt.figure(figsize=(args.figure_width, args.figure_height))
    render_all(
        fig,
        out_dir,
        method_order,
        color_palette,
        read_depth=args.read_depth
    )

    draw_legend(fig, method_order, color_palette)

    out_path = out_dir / f'errors_{read_depth}.{args.format}'

    plt.savefig(out_path)
    fig.tight_layout()
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    main()
