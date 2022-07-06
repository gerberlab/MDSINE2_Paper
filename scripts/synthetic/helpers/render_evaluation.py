from pathlib import Path
import argparse

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help='<Required> The target output directory. Should contain the results of `evaluate` script.')
    parser.add_argument('-f', '--format', type=str, required=False,
                        default='pdf',
                        help='<Optional> The plot image format (Default: pdf)')
    return parser.parse_args()


def load_df(df_path: Path) -> pd.DataFrame:
    df = pd.read_csv(df_path)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    df = df.sort_values(['ReadDepth', 'NoiseLevel'], ascending=[True, True])
    return df


def render_growth_rate_errors(df: pd.DataFrame, text_ax, ax1, ax2, order, palette):
    df['x'] = df['ReadDepth'].astype(str) + ' reads\n' + df['NoiseLevel'].astype(str) + ' noise'
    sb.barplot(
        data=df.loc[df['ReadDepth'] == 1000],
        x='x',
        y='Error',
        hue='Method',
        ax=ax1,
        hue_order=order,
        palette=palette
    )
    sb.barplot(
        data=df.loc[df['ReadDepth'] == 25000],
        x='x',
        y='Error',
        hue='Method',
        ax=ax2,
        hue_order=order,
        palette=palette
    )
    ax1.set_ylabel('RMSE')
    ax2.set_ylabel(None)
    ax1.set_xlabel(None)
    ax2.set_xlabel(None)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    text_ax.text(x=0.5, y=0.5, s='Growth Rates', fontsize=14, horizontalalignment='center', verticalalignment='center')
    text_ax.axis('off')


def render_interaction_strength_errors(df: pd.DataFrame, text_ax, ax1, ax2, order, palette):
    df['x'] = df['ReadDepth'].astype(str) + ' reads\n' + df['NoiseLevel'].astype(str) + ' noise'
    sb.barplot(
        data=df.loc[df['ReadDepth'] == 1000],
        x='x',
        y='Error',
        hue='Method',
        ax=ax1,
        hue_order=order,
        palette=palette
    )
    sb.barplot(
        data=df.loc[df['ReadDepth'] == 25000],
        x='x',
        y='Error',
        hue='Method',
        ax=ax2,
        hue_order=order,
        palette=palette
    )
    ax1.set_ylabel('RMSE')
    ax2.set_ylabel(None)
    ax1.set_xlabel(None)
    ax2.set_xlabel(None)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    text_ax.text(x=0.5, y=0.5, s='Interaction Strengths', fontsize=14, horizontalalignment='center', verticalalignment='center')
    text_ax.axis('off')


def render_holdout_trajectory_errors(df: pd.DataFrame, text_ax, ax1, ax2, order, palette):
    df['x'] = df['ReadDepth'].astype(str) + ' reads\n' + df['NoiseLevel'].astype(str) + ' noise'
    sb.barplot(
        data=df.loc[df['ReadDepth'] == 1000],
        x='x',
        y='Error',
        hue='Method',
        ax=ax1,
        hue_order=order,
        palette=palette
    )
    sb.barplot(
        data=df.loc[df['ReadDepth'] == 25000],
        x='x',
        y='Error',
        hue='Method',
        ax=ax2,
        hue_order=order,
        palette=palette
    )
    ax1.set_ylabel('RMSE')
    ax2.set_ylabel(None)
    ax1.set_xlabel(None)
    ax2.set_xlabel(None)
    text_ax.text(x=0.5, y=0.5, s='Holdout trajectory', fontsize=14, horizontalalignment='center', verticalalignment='center')
    text_ax.axis('off')


def render_topology_errors(df: pd.DataFrame, text_ax, ax1, ax2, order, palette):
    def auroc(_df):
        _df = _df.sort_values('FPR', ascending=True)
        fpr = np.concatenate([[0.], _df['FPR'].to_numpy(), [1.]])
        tpr = np.concatenate([[0.], _df['TPR'].to_numpy(), [1.]])
        return scipy.integrate.trapz(
            y=tpr,
            x=fpr,
        )

    area_df = df.groupby(['Method', 'ReadDepth', 'NoiseLevel']).apply(auroc).rename('Error').reset_index()

    area_df['x'] = area_df['ReadDepth'].astype(str) + ' reads\n' + area_df['NoiseLevel'].astype(str) + ' noise'
    sb.barplot(
        data=area_df.loc[area_df['ReadDepth'] == 1000],
        x='x',
        y='Error',
        hue='Method',
        ax=ax1,
        hue_order=order,
        palette=palette
    )
    sb.barplot(
        data=area_df.loc[area_df['ReadDepth'] == 25000],
        x='x',
        y='Error',
        hue='Method',
        ax=ax2,
        hue_order=order,
        palette=palette
    )

    ax1.set_ylabel('AUROC')
    ax2.set_ylabel(None)
    ax1.set_xlabel(None)
    ax2.set_xlabel(None)
    text_ax.text(x=0.5, y=0.5, s='Network Structure', fontsize=14, horizontalalignment='center', verticalalignment='center')
    text_ax.axis('off')


def render_all(dataframe_dir: Path, output_path: Path):
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    method_order = ['MDSINE2', 'MDSINE1', 'lra-elastic_net', 'glv-elastic_net', 'glv-ridge', 'glv-ra-elastic_net', 'glv-ra-ridge']
    palette = {method: c for method, c in zip(method_order, default_colors)}

    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    spec = fig.add_gridspec(ncols=4, nrows=4, height_ratios=[1, 10, 1, 10], width_ratios=[1, 1, 1, 1])

    ax0 = fig.add_subplot(spec[0, :2])
    ax1, ax2 = fig.add_subplot(spec[1, 0]), fig.add_subplot(spec[1, 1])
    render_growth_rate_errors(
        load_df(dataframe_dir / "growth_rate_errors.csv"),
        text_ax=ax0,
        ax1=ax1,
        ax2=ax2,
        order=method_order,
        palette=palette
    )

    ax0 = fig.add_subplot(spec[0, 2:])
    ax1, ax2 = fig.add_subplot(spec[1, 2]), fig.add_subplot(spec[1, 3])
    render_interaction_strength_errors(
        load_df(dataframe_dir / "interaction_strength_errors.csv"),
        text_ax=ax0,
        ax1=ax1,
        ax2=ax2,
        order=method_order,
        palette=palette
    )

    ax0 = fig.add_subplot(spec[2, :2])
    ax1, ax2 = fig.add_subplot(spec[3, 0]), fig.add_subplot(spec[3, 1])
    render_holdout_trajectory_errors(
        load_df(dataframe_dir / 'holdout_trajectory_errors.csv'),
        text_ax=ax0,
        ax1=ax1,
        ax2=ax2,
        order=method_order,
        palette=palette
    )

    ax0 = fig.add_subplot(spec[2, 2:])
    ax1, ax2 = fig.add_subplot(spec[3, 2]), fig.add_subplot(spec[3, 3])
    render_topology_errors(
        load_df(dataframe_dir / "topology_errors.csv"),
        text_ax=ax0,
        ax1=ax1,
        ax2=ax2,
        order=method_order,
        palette=palette
    )

    plt.savefig(output_path)


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    render_all(out_dir, out_dir / f'errors.{args.format}')


if __name__ == "__main__":
    main()
