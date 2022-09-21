import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import argparse
import pandas as pd
import mdsine2 as md2
from mdsine2.names import STRNAMES

from matplotlib.colors import LogNorm
from pathlib import Path


def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument("--base_loc", help="directory where necessary files to make"
                                           "the plot are saved")
    parser.add_argument("--module_result_pkl", help="name of the pkl file "
                        "containing the results of enrichment analysis on modules")
    parser.add_argument("--pathways_result_pkl", help="name of the pkl file "
                        "containing the results of enrichment analysis on pathways")
    parser.add_argument("--cazyme_result_pkl", help="name of the pkl file"
                        "containing the results of enrichment analysis on cazymes")
    parser.add_argument("--output_loc", help="location where the output is saved")
    parser.add_argument("--mcmc_file_loc", help="location of the mcmc pkl files")
    return parser.parse_args()


def load_pkl_file(file_loc):
    """
    :param (Path) file_loc: loads the pkl file and returns the data saved in the file
    """
    print(file_loc)
    file = pickle.load(open(file_loc, "rb"))

    return file


def plot_p_value_heatmap(axes, df, title, y_label, plot_cbar=False, cbar_ax=None,
                         vmax=0.5, vmin=0.0005, label_x=False):
    """
    plots the p values as a heatmpa
    """
    cmap = sns.color_palette("Blues", n_colors=5)

    np_df = df.to_numpy()
    np_df = np.where(np_df>0.05, 1.5, np_df)
    np_df = np.where((0.01<np_df) & (np_df <=0.05), 2.5, np_df)
    np_df = np.where((0.001<np_df) & (np_df <=0.01), 3.5, np_df)
    np_df = np.where((0.0001<np_df) & (np_df <=0.001), 4.5, np_df)
    np_df = np.where(np_df <=0.0001, 5.5, np_df)
    df_new = pd.DataFrame(np_df, columns=df.columns, index=df.index)

    if plot_cbar:
        heatmap = sns.heatmap(df_new, ax=axes, cmap=cmap, linewidths=0.5,
                              yticklabels=True, xticklabels=True,
                              cbar_ax=cbar_ax, vmax=6, vmin=1,
                              cbar_kws={"ticks":[1.5, 2.5, 3.5, 4.5, 5.5]})
        cbar = heatmap.collections[0].colorbar
        cbar.set_ticklabels(["ns", "$p\\leq0.05$", "$p\\leq0.01$", "$p\\leq0.001$",
                             "$p\\leq0.0001$"])
        cbar.ax.tick_params(labelsize=40, labelleft =False, labelright=True,
                            left=False, right=True, length=20, width=1)
        cbar.ax.tick_params(length=0, width=0, which="minor")
        cbar.ax.set_title("Adjusted\n p-values\n", size=40,
                          fontweight="bold")
    else:
        sns.heatmap(df_new, ax=axes, cmap=cmap, linewidths=0.5, yticklabels=True,
                    cbar=False, xticklabels=True, vmax=6, vmin=1)
    axes.set_xlabel("Module", fontsize=25, fontweight="bold", labelpad=2)
    axes.set_ylabel(y_label, fontsize=20, fontweight="bold", labelpad=10)
    axes.set_title(title, loc="left", fontsize=30, fontweight="bold")
    axes.set_yticklabels(axes.get_yticklabels(), rotation=0, fontsize=20)
    axes.tick_params(length=0, axis="both")
    xtick_labels = []
    cohort_label="M"
    xtick_labels = [cohort_label + str(i) for i in range(1,
        df.to_numpy().shape[1] + 1)]
    axes.set_xticklabels(xtick_labels, fontsize=20, rotation=90)
    axes.tick_params(length=0, axis="both")
    axes.set_aspect('equal')

    return axes


def format_df(df, n_cluster):
    columns = df.columns
    shape = df.shape

    missing = [1] * shape[0]
    for i in range(1, n_cluster + 1):
        if "Cluster " + str(i) not in columns:
            df["Cluster " + str(i)] = missing
    sorted_cols = ["Cluster " + str(i) for i in range(1, n_cluster + 1)]
    new_df = df[sorted_cols]
    cols = [str(i) for i in range(1, n_cluster + 1)]
    new_df.columns = cols

    return new_df


if __name__ == "__main__":

    args = parse_arguments()
    base_loc = Path(args.base_loc)
    mcmc = md2.BaseMCMC.load(args.mcmc_file_loc)
    clustering = mcmc.graph[STRNAMES.CLUSTERING].clustering
    cluster_size = len(clustering)
    cazyme_results = format_df(load_pkl_file(base_loc/"{}.pkl".format(
        args.cazyme_result_pkl)), cluster_size)
    modules_results = format_df(load_pkl_file(base_loc/"{}.pkl".format(
        args.module_result_pkl)), cluster_size)
    pathways_results = format_df(load_pkl_file(base_loc/"{}.pkl".format(
        args.pathways_result_pkl)), cluster_size)

    fig = plt.figure(figsize=(20, 12))
    spec = gridspec.GridSpec(ncols=8, nrows=11, figure=fig)
    axes1 = fig.add_subplot(spec[0:5, 1:7])
    axes2 = fig.add_subplot(spec[5:8, 2:8])
    axes3 = fig.add_subplot(spec[8:11, 2:8])

    cbar_axes = fig.add_subplot(spec[4:9, 7])
    plot_p_value_heatmap(axes1, modules_results, "KEGG Modules", "Module")
    plot_p_value_heatmap(axes2, pathways_results, "KEGG Pathways", "Module",
                         plot_cbar=True, cbar_ax=cbar_axes)
    plot_p_value_heatmap(axes3, cazyme_results, "CaZyme", "Module")

    output_loc = Path(args.output_loc)
    output_loc.mkdir(parents=True, exist_ok=True)

    plt.savefig(output_loc/"functional_enrichment.pdf", bbox_inches="tight")