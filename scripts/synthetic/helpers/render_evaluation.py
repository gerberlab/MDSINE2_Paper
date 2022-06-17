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


def load_df(df_path: Path) -> pd.DataFrame:
    df = pd.read_csv(df_path)
    df['NoiseLevel'] = pd.Categorical(
        df['NoiseLevel'],
        categories=['low', 'medium', 'high']
    )
    df = df.sort_values(['ReadDepth', 'NoiseLevel'], ascending=[True, True])
    return df


def render_growth_rate_errors(df: pd.DataFrame, ax):
    df['x'] = df['ReadDepth'].astype(str) + ' reads\n' + df['NoiseLevel'].astype(str) + ' noise'
    sb.barplot(
        data=df,
        x='x',
        y='Error',
        hue='Method',
        ax=ax
    )
    # sb.swarmplot(
    #     data=df,
    #     x='x',
    #     y='Error',
    #     hue='Method',
    #     dodge=True,
    #     ax=ax
    # )
    ax.set_ylabel('RMSE')


def render_interaction_strength_errors(df: pd.DataFrame, ax):
    df['x'] = df['ReadDepth'].astype(str) + ' reads\n' + df['NoiseLevel'].astype(str) + ' noise'
    sb.barplot(
        data=df,
        x='x',
        y='Error',
        hue='Method',
        ax=ax
    )
    # sb.swarmplot(
    #     data=df,
    #     x='x',
    #     y='Error',
    #     hue='Method',
    #     dodge=True,
    #     ax=ax
    # )
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

    df['x'] = df['ReadDepth'].astype(str) + ' reads\n' + df['NoiseLevel'].astype(str) + ' noise'
    area_df = df.groupby(['Method', 'ReadDepth', 'NoiseLevel']).apply(auroc)
    print(area_df)


def render_all(dataframe_dir: Path, output_path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    render_growth_rate_errors(
        load_df(dataframe_dir / "growth_rate_errors.csv"),
        axes[0, 0]
    )
    render_interaction_strength_errors(
        load_df(dataframe_dir / "interaction_strength_errors.csv"),
        axes[1, 0]
    )
    render_topology_errors(
        load_df(dataframe_dir / "topology_errors.csv"),
        axes[1, 1]
    )

    plt.savefig(output_path)


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    render_all(out_dir, out_dir / f'errors.{args.format}')


if __name__ == "__main__":
    main()
