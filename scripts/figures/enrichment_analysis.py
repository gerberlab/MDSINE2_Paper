import numpy as np
import pandas as pd
import copy
from scipy.special import comb
from scipy.stats import hypergeom
from statsmodels.stats import multitest as mtest
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import os
import pickle
import mdsine2 as md2
import math
from mdsine2.names import STRNAMES
from matplotlib.gridspec import GridSpec
from pandas.plotting import table
from matplotlib.colors import ListedColormap
import argparse
#import pylab as pl

from matplotlib import rcParams
from matplotlib import font_manager

rcParams['pdf.fonttype'] = 42

font_dirs = ['gibson_inference/figures/arial_fonts']
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
#print(font_files)

for font_file in font_files:
    #print(font_file)
    ff = font_file.split("/")[-1]
    if "._" not in ff:
        font_manager.fontManager.addfont(font_file)

# change font
rcParams['font.family'] = 'Arial'

def compute_p_value(N, M, n, k):
    '''computes and returns the hypergeometric p-values
       @Parameters
       ------------------------------------------------------------------
       N, M, n, k: (int)

       @returns
        ------------------------------------------------------------------------
       float
    '''
    lim = min(n, M)
    p = 0
    for i in range(k, lim + 1):
        p += comb(M, i) * comb(N - M, n - i) / comb(N, n)

    p1 = hypergeom.sf(k-1, N, M, n) #to check the accuracy of p

    return p

def obtain_p_vals(cluster_dict, hierarchy_otu_dict):
    total_otus = sum([len(cluster_dict[id_]) for id_ in cluster_dict])

    cluster_all_p ={}
    cluster_all_level = {}
    for id_ in sorted(cluster_dict.keys()):
        #print("id:", id_)
        otus_li = cluster_dict[id_]
        cluster_all_p[id_] = []
        cluster_all_level[id_] = []
        n_otus_cluster = len(otus_li)

        for level in hierarchy_otu_dict:
            level_li = hierarchy_otu_dict[level]
            n_otus_annotated = len(level_li)
            n_otus_annotated_cluster = len(set(otus_li).intersection(level_li))
            if n_otus_annotated_cluster != 0 :
                p_val = compute_p_value(total_otus, n_otus_annotated,
                    n_otus_cluster, n_otus_annotated_cluster)
                cluster_all_p[id_].append(p_val)
                cluster_all_level[id_].append(level)

    return cluster_all_p, cluster_all_level

def clusterize(labels, taxa_list):
    """returns the consenus cluster as a dictionary ((int) cluster_id ->
    ([str]) ids of OTUs belonging to the cluster)

    @parameters
    labels : ([int]) the consensus cluster id of OTUs in taxa_list
    taxa_list : ([str]) ids of OTU
    """

    cluster = {}
    for i in range(len(labels)):
        if labels[i] + 1 not in cluster:
            cluster[labels[i] + 1] = []
        cluster[labels[i] + 1].append(taxa_list[i])

    return cluster

def map_level_to_otu(taxas_info, level):
    map_ = {}
    for taxa in taxas_info:
        key = taxa.taxonomy[level]
        if key not in map_:
            map_[key] = []
        name = taxa.name
        map_[key].append(name)

    return map_

def parse_cluster(mcmc):

    clustering = mcmc.graph[STRNAMES.CLUSTERING].clustering
    consensus_cluster_labels = clustering.toarray()
    taxa_list = []
    taxas = mcmc.graph.data.taxa
    for taxa in taxas:
        taxa_list.append(taxa.name)

    cluster = clusterize(consensus_cluster_labels, taxa_list)
    #export_cluster(cluster, taxa_list, study_obj.taxa, name)
    return cluster

def pivot_df(p_li, p_info, adjusted_p_li, adjusted_p_res_li, name):

    cluster_row = []
    hierarchy_names = []
    values = []
    enriched = []

    for i in range(len(p_li)):
        if adjusted_p_res_li[i] == True:
            cluster_row.append(" ".join(p_info[i].split()[:2]))
            hierarchy_names.append(p_info[i].split()[-1])
            values.append(adjusted_p_li[i])
            enriched.append(adjusted_p_res_li[i])

    data_frame = pd.DataFrame({"cluster_id":cluster_row,
            "module_name":hierarchy_names, "p_value":values})
    df_pivot = data_frame.pivot(values="p_value", index="module_name",
                columns="cluster_id")
    df_pivot = df_pivot.fillna(1)

    #df_pivot.to_csv(name + ".csv", sep=",")
    return df_pivot

def run_enrichment(mcmc, hierarchy_level, pivot_name):

    cluster_dict = parse_cluster(mcmc)
    taxas = mcmc.graph.data.taxa
    level_otu_map = map_level_to_otu(taxas, hierarchy_level)

    p_values, p_module_info = obtain_p_vals(cluster_dict, level_otu_map)

    all_p_li = []
    all_p_info_li = []

    for keys in sorted(p_values.keys()):
        all_p_li = all_p_li + p_values[keys]
        modules_cluster = []
        for key2 in p_module_info[keys]:
            modules_cluster.append("Cluster " + str(keys) + " " + key2)
        all_p_info_li = all_p_info_li + modules_cluster

    adjusted_p = []

    if len(all_p_li) != 0:
        adjusted_p = mtest.multipletests(copy.deepcopy(all_p_li), alpha=0.05,
            method="fdr_bh", is_sorted=False)
    else:
        print("There are no valid p values ")

    pivoted_df = pivot_df(all_p_li, all_p_info_li, adjusted_p[1], adjusted_p[0],
        pivot_name)

    df_copy = copy.deepcopy(pivoted_df)
    index = np.sort([int(idx.split()[1]) for idx in df_copy.columns])
    index = ["Cluster " + str(i) for i in index]
    df_save = df_copy[index]
    #export_df(cluster_dict, level_otu_map, all_p_li, all_p_info_li,
        #adjusted_p[1], pivot_name)

    #loc = "gibson_inference/figures/output_figures/p_supplemental_figure5_"
    #os.makedirs("{}".format(loc), exist_ok=True)
    #df_save.to_csv("{}/p_test_{}.csv".format(loc, pivot_name), sep=",")

    return pivoted_df

def export_df(cluster_info, level_otu_dict, raw_p, p_info, adj_p, savename):

    array_adj_p = np.ones((len(level_otu_dict), len(cluster_info)))
    array_p = np.ones((len(level_otu_dict), len(cluster_info)))
    index_dict = {}
    i = 0
    for key in np.sort(list(level_otu_dict.keys())):
        index_dict[key] = i
        i += 1

    N = len(raw_p)
    print("check:", len(raw_p), len(p_info), len(adj_p))
    for n in range(N):
        info_row = p_info[n].split()
        j = int(info_row[1]) - 1
        i = index_dict[info_row[2]]
        array_p[i, j] = raw_p[n]
        array_adj_p[i, j] = adj_p[n]
    loc = "gibson_inference/figures/output_figures/p_supplemental_figure5"
    os.makedirs("{}".format(loc), exist_ok=True)

    combined_array = np.ones((len(level_otu_dict), len(cluster_info)*2))
    j = 0
    for i in range(len(cluster_info)):
        combined_array[:, j] = array_p[:, i]
        combined_array[:, j+1] = array_adj_p[:, i]
        j += 2

    df_p = pd.DataFrame(array_p)
    df_p.index = np.sort(list(level_otu_dict.keys()))
    df_p.columns = ["Cluster " + str(i+1) for i in range(len(cluster_info))]

    df_adj_p = pd.DataFrame(array_adj_p)
    df_adj_p.index = np.sort(list(level_otu_dict.keys()))
    df_adj_p.columns = ["Cluster " + str(i+1) for i in range(len(cluster_info))]

    df_combined = pd.DataFrame(combined_array)
    df_combined.index = np.sort(list(level_otu_dict.keys()))
    df_combined.columns = ["Raw p-value", "BH adjusted p-value"] * len(cluster_info)

    df_p.to_csv("{}/p_{}.csv".format(loc, savename), sep=",")
    df_adj_p.to_csv("{}/adj_p_{}.csv".format(loc, savename), sep=",")
    df_combined.to_csv("{}/combined_p_{}.csv".format(loc, savename), sep=",")


def pivot_cluster_membership(mcmc, level):

    taxas = mcmc.graph.data.taxa
    level_otu_map = map_level_to_otu(taxas, level)
    cluster_dict = parse_cluster(mcmc)

    data = []
    index = []
    for level in sorted(level_otu_map.keys()):
        name_li = level.split("_")
        if len(name_li) != 1:
            name_li = name_li[:-1]
        new_name = " ".join(name_li)

        index.append(new_name)
        level_info_row = []
        for id_ in sorted(cluster_dict.keys()):
            level_info_row.append(len(set(cluster_dict[id_]).intersection(
                level_otu_map[level])))
        data.append(level_info_row)

    columns = ["{}".format(i) for i in range(1, len(cluster_dict) +1)]
    df = pd.DataFrame(np.asarray(data), index=index, columns=columns)

    return df

def plot_module(cohort, df, axes, title, ylab, plot_cbar=False, cbar_ax=None, vmax=0.5,
    vmin=0.0005, label_x=False, label_x_gram=False, labels_gram_stain=None,
    tick_right=False):

    cmap = sns.color_palette("Blues", n_colors=5)

    np_df = df.to_numpy()
    np_df = np.where(np_df>0.05, 1.5, np_df)
    np_df = np.where((0.01<np_df) & (np_df <=0.05), 2.5, np_df)
    np_df = np.where((0.001<np_df) & (np_df <=0.01), 3.5, np_df)
    np_df = np.where((0.0001<np_df) & (np_df <=0.001), 4.5, np_df)
    np_df = np.where(np_df <=0.0001, 5.5, np_df)

    df_new = pd.DataFrame(np_df, columns=df.columns, index=df.index)


    if plot_cbar:
        map_ = sns.heatmap(df_new, ax=axes, cmap=cmap,
        linewidths=0.5, yticklabels=True, xticklabels=True,
        cbar_ax=cbar_ax,vmax=6, vmin=1, cbar_kws={"ticks":[1.5, 2.5, 3.5,
        4.5, 5.5]})
        cbar = map_.collections[0].colorbar
        #cbar.set_ticklabels(["p<0.0005", "p<0.005", "p<0.05", "p>0.05"])
        cbar.set_ticklabels(["ns",  "$p\\leq0.05$", "$p\\leq0.01$", "$p\\leq0.001$",
            "$p\\leq0.0001$"])
        cbar = map_.collections[0].colorbar
        cbar.ax.tick_params(labelsize=40, labelleft =False, labelright = True,
            left = False, right = True, length = 20, width = 1)
        cbar.ax.tick_params(length = 0, width = 0, which = "minor")
        cbar.ax.set_title("Adjusted\n p-values\n", size = 40,
             fontweight = "bold")

    else:
        sns.heatmap(df_new, ax=axes, cmap=cmap, linewidths=0.5, yticklabels=True,
         cbar=False, xticklabels=True,vmax=6, vmin=1)


    axes.set_xlabel("Module", fontsize=45, fontweight="bold", labelpad=5)
    axes.set_ylabel(ylab, fontsize=35, fontweight="bold", labelpad=25)
    axes.set_title(title, loc="left", fontsize=50, fontweight="bold")
    axes.set_yticklabels(axes.get_yticklabels(), rotation = 0,
               fontsize = 35)
    axes.tick_params(length=0, axis = "both")
    xtick_labels = []
    cohort_label = "M"
    if cohort != "healthy":
        cohort_label = "M"
    xtick_labels = [cohort_label + str(i) for i in range(1,
        df.to_numpy().shape[1] + 1)]

    axes.set_xticklabels(xtick_labels, fontsize = 35, rotation=90)
    #if label_x:
    #    axes.set_xticklabels(axes.get_xticklabels(), rotation = 90, fontsize = 50)
    #else:
    #    axes.set_xticklabels([], rotation = 0, fontsize = 18)

    if tick_right:
        axes.tick_params(length=0, axis = "both", labelright=True, labelleft=False)
        axes.yaxis.set_label_position("right")
    else:
        axes.tick_params(length=0, axis = "both")

    axes.set_aspect('equal')
    return axes

def _make_table(cohort, table_df, axes, title, tick_right, y_lab):

    map_ = sns.heatmap(table_df, linewidths=1, ax=axes, cmap=ListedColormap(['white']),
        cbar=False, annot=True, linecolor="black",
        annot_kws={"fontsize":30})
    cohort_label = "M"
    if cohort != "healthy":
        cohort_label = "M"
    xtick_labels = [cohort_label + str(i) for i in range(1,
        table_df.to_numpy().shape[1] + 1)]

    #axes.set_xlabel(x_lab, fontsize=35)
    axes.set_yticklabels(axes.get_yticklabels(), rotation = 0,
               fontsize = 35)
    axes.set_xticklabels(xtick_labels, rotation = 90,
               fontsize = 35)
    axes.set_ylabel(y_lab, fontsize=35, fontweight="bold", labelpad=25)
    axes.set_xlabel("Module", fontsize=45, fontweight="bold", labelpad=5)
    axes.set_title(title, loc="left", fontsize=50, fontweight="bold")

    if tick_right:
        axes.tick_params(length=0, axis = "both", labelright=True, labelleft=False)
        axes.yaxis.set_label_position("right")
    else:
        axes.tick_params(length=0, axis = "both")
    return axes

def make_plot(df_healthy_p_order, df_healthy_p_family, df_healthy_p_class,
              df_healthy_p_phylum, df_healthy_order_abund,
              df_healthy_family_abund,  df_healthy_class_abund,
              df_healthy_phylum_abund, loc):

    fig = plt.figure(figsize=(40, 70))
    spec = GridSpec(ncols=44, nrows= 110, figure=fig)

    family_abund_start = 2
    family_abund_end = 26

    order_abund_start = 29
    order_abund_end = 41

    class_abund_start = 44
    class_abund_end = 56

    phylum_abund_start = 59
    phylum_abund_end = 65

    family_enrich_start = 70
    family_enrich_end = 77

    order_enrich_start = 80
    order_enrich_end = 88

    class_enrich_start = 91
    class_enrich_end = 99

    phylum_enrich_start = 102
    phylum_enrich_end = 110

    healthy_column_start = 2
    healthy_column_end = 18

    cax_row_start = 78
    cax_row_end = 92

    cax_column_start = 21
    cax_column_end = 22

    ax_healthy_family_abund = fig.add_subplot(spec[family_abund_start:family_abund_end,
       healthy_column_start:healthy_column_end])

    ax_healthy_order_abund = fig.add_subplot(spec[order_abund_start:order_abund_end,
       healthy_column_start:healthy_column_end])

    ax_healthy_class_abund = fig.add_subplot(spec[class_abund_start:class_abund_end,
       healthy_column_start:healthy_column_end])

    ax_healthy_phylum_abund = fig.add_subplot(spec[phylum_abund_start:phylum_abund_end,
       healthy_column_start:healthy_column_end])

    ax_healthy_family_enrich = fig.add_subplot(spec[family_enrich_start:family_enrich_end,
       healthy_column_start:healthy_column_end])

    ax_healthy_order_enrich = fig.add_subplot(spec[order_enrich_start:order_enrich_end,
       healthy_column_start:healthy_column_end])

    ax_healthy_class_enrich = fig.add_subplot(spec[class_enrich_start:class_enrich_end,
       healthy_column_start:healthy_column_end])

    ax_healthy_phylum_enrich = fig.add_subplot(spec[phylum_enrich_start:phylum_enrich_end,
       healthy_column_start:healthy_column_end])


    cax = fig.add_subplot(spec[cax_row_start:cax_row_end,
        cax_column_start:cax_column_end])

    ax_healthy_family_enrich = plot_module("healthy", df_healthy_p_family, ax_healthy_family_enrich,
        "E", "Family", plot_cbar=False)
    ax_healthy_order_enrich = plot_module("healthy",df_healthy_p_order, ax_healthy_order_enrich,
        "F", "Order", plot_cbar=False)
    ax_healthy_class_enrich = plot_module("healthy", df_healthy_p_class, ax_healthy_class_enrich,
        "G", "Class", plot_cbar=True, cbar_ax=cax)
    ax_healthy_phylum_enrich = plot_module("healthy", df_healthy_p_phylum, ax_healthy_phylum_enrich,
        "H", "Phylum", plot_cbar=False)
    ax_healthy_family_abund = _make_table("healthy", df_healthy_family_abund,
        ax_healthy_family_abund, "A", False, "Family")
    ax_healthy_order_abund = _make_table("healthy", df_healthy_order_abund,
        ax_healthy_order_abund, "B", False, "Order")
    ax_healthy_class_abund = _make_table("healthy", df_healthy_class_abund,
        ax_healthy_class_abund, "C", False, "Class")
    ax_healthy_phylum_abund = _make_table("healthy", df_healthy_phylum_abund,
        ax_healthy_phylum_abund, "D", False, "Phylum")

    fig.text(0.025, 0.075, "NA: Taxonomy not resolved", fontsize=45)
    fig.subplots_adjust(hspace=25)

    os.makedirs(loc, exist_ok=True)
    fig.text(0.255, 0.88, "Healthy", fontsize=55, fontweight="bold")
    fig.savefig(loc + "/enrichment.pdf", bbox_inches="tight")


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

def parse_args():

    parser = argparse.ArgumentParser(description = "files needed for making"\
    "supplemental figure 5")
    parser.add_argument("-loc1", "--healthy_mcmc_loc", required = "True",
       help = "a pl.BaseMCMC pkl file for healthy runs")
    parser.add_argument("-o_loc", "--output_loc", required="True",
        help = "directory(folder name) where the output figure is saved")

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()
    print("Making Supplemental Figure 5")

    mcmc_healthy = md2.BaseMCMC.load(args.healthy_mcmc_loc)

    healthy_order_enrichment = run_enrichment(mcmc_healthy, "order", "healthy_order")
    healthy_family_enrichment = run_enrichment(mcmc_healthy, "family", "healthy_family")
    healthy_class_enrichment = run_enrichment(mcmc_healthy, "class", "healthy_class")
    healthy_phylum_enrichment = run_enrichment(mcmc_healthy, "phylum", "healthy_phylum")
    healthy_family_abundance = pivot_cluster_membership(mcmc_healthy, "family")
    healthy_order_abundance = pivot_cluster_membership(mcmc_healthy, "order")
    healthy_class_abundance = pivot_cluster_membership(mcmc_healthy, "class")
    healthy_phylum_abundance = pivot_cluster_membership(mcmc_healthy, "phylum")

    healthy_N = healthy_order_abundance.to_numpy().shape[1]

    healthy_order_enrichment = format_df(healthy_order_enrichment, healthy_N)
    healthy_family_enrichment= format_df(healthy_family_enrichment, healthy_N)
    healthy_class_enrichment= format_df(healthy_class_enrichment, healthy_N)
    healthy_phylum_enrichment= format_df(healthy_phylum_enrichment, healthy_N)

    make_plot(healthy_order_enrichment, healthy_family_enrichment, healthy_class_enrichment,
             healthy_phylum_enrichment, healthy_order_abundance,  healthy_family_abundance,
             healthy_class_abundance, healthy_phylum_abundance,  args.output_loc)

    print("Done Making Supplemental Figure 5")
