from pathlib import Path
import argparse

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


def render_growth_rate_errors(df: pd.DataFrame, ax):
    df = df.sort_values(['ReadDepth', 'NoiseLevel'], ascending=[True, True])
    df['x'] = df['ReadDepth'].astype(str) + ' reads, ' + df['NoiseLevel'].astype(str) + ' noise'
    sb.boxplot(
        data=df,
        x='x',
        y='Error',
        hue='Method',
        showfliers=False,
        ax=ax
    )
    sb.swarmplot(
        data=df,
        x='x',
        y='Error',
        hue='Method',
        dodge=True,
        ax=ax
    )
    ax.set_ylabel('RMSE')


def render_interaction_strength_errors(df: pd.DataFrame, ax):
    df = df.sort_values(['ReadDepth', 'NoiseLevel'], ascending=[True, True])
    df['x'] = df['ReadDepth'].astype(str) + ' reads, ' + df['NoiseLevel'].astype(str) + ' noise'
    sb.boxplot(
        data=df,
        x='x',
        y='Error',
        hue='Method',
        showfliers=False,
        ax=ax
    )
    sb.swarmplot(
        data=df,
        x='x',
        y='Error',
        hue='Method',
        dodge=True,
        ax=ax
    )
    ax.set_ylabel('RMSE')


def render_topology_errors(df: pd.DataFrame, ax):
    def auroc(_df):
        _df = _df.sort_values('FPR', ascending=True)
        fpr = _df['FPR']
        tpr = _df['TPR']
        return scipy.integrate.trapz(
            y=tpr,
            x=fpr,
        )

    df = df.sort_values(['ReadDepth', 'NoiseLevel'], ascending=[True, True])
    df['x'] = df['ReadDepth'].astype(str) + ' reads, ' + df['NoiseLevel'].astype(str) + ' noise'
    area_df = df.groupby(['Method', 'ReadDepth', 'NoiseLevel']).apply(auroc)
    print(area_df)


def render_all(dataframe_dir: Path, output_path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))

    render_growth_rate_errors(
        pd.read_csv(dataframe_dir / "growth_rate_errors.csv"),
        axes[0, 0]
    )
    render_interaction_strength_errors(
        pd.read_csv(dataframe_dir / "interaction_strength_errors.csv"),
        axes[1, 0]
    )
    render_topology_errors(
        pd.read_csv(dataframe_dir / "topology_errors"),
        axes[1, 1]
    )

    plt.savefig(output_path)


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    render_all(out_dir, out_dir / f'errors.{args.format}')


if __name__ == "__main__":
    main()
