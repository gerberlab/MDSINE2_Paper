#figure that plots both box plots and abundance plot in the same figure
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats
import pickle as pkl
import pandas as pd
import seaborn as sns
import argparse
import matplotlib.gridspec as gridspec
import mdsine2 as md2

import matplotlib.lines as mlines
from statsmodels.stats import multitest
from matplotlib import rcParams
from matplotlib import font_manager
from pathlib import Path

TITLE_FONTSIZE = 20
TICK_FONTSIZE = 16
AXES_FONTSIZE = 18

TITLE_FONTSIZE = 20
TICK_FONTSIZE = 16
AXES_FONTSIZE = 18

REL_ORDER = ["MDSINE2", "cLV", "LRA", "gLV-RA", "gLV-ridge", "gLV-elastic\n net"]
ABS_ORDER = ["MDSINE2", "gLV-ridge", "gLV-elastic\n net"]

HEX_REL = sns.color_palette("tab10").as_hex()
HEX_ABS = sns.color_palette("tab10").as_hex()

PAL_REL = {"MDSINE2":HEX_REL[0], "cLV":HEX_REL[3], "LRA":HEX_REL[4],
           "gLV-RA":HEX_REL[5], "gLV-ridge":HEX_REL[1],
           "gLV-elastic\n net":HEX_REL[2]}
PAL_ABS = {"MDSINE2":HEX_REL[0], "gLV-ridge":HEX_REL[1],
           "gLV-elastic\n net":HEX_REL[2]}

rcParams['pdf.fonttype'] = 42
font_dirs = ['figures/arial_fonts']
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

for font_file in font_files:
    ff = font_file.split("/")[-1]
    if "._" not in ff:
        font_manager.fontManager.addfont(font_file)
# change font
rcParams['font.family'] = 'Arial'


def parse_args():

    parser = argparse.ArgumentParser(description = "files needed for making"\
                                                   "figure 3")
    parser.add_argument("-path1", "--mdsine_path", required = True,
                        help = "path lo mdsine2 forward sims")
    parser.add_argument("-path2", "--non_glv_elas_path", required = True,
                        help = "path to all model sims but gLV trained using elastic net")
    parser.add_argument("-path3", "--non_glv_ridge_path", required = True,
                        help = "path to all models but gLV sims trained using ridge net")
    parser.add_argument("-path4", "--glv_elas_path", required = True,
                        help = "path to gLV model sims trained using elastic net")
    parser.add_argument("-path5", "--glv_ridge_path", required = True,
                        help = "path to gLV model sims trained using ridge")
    parser.add_argument("-o_path", "--output_path", required=True,
                        help = "path to the folder where the output is saved")
    parser.add_argument("-o_name", "--output_name", required=True)
    parser.add_argument("-sf", "--study_file", required=True)
    parser.add_argument("-rl", "--rel_abund_lim", required=True, type=float)

    return parser.parse_args()


def compute_rms_error(pred, truth):
    """
    computes the root mean square error
    (np.ndarray) pred, truth: arrays containing true and predicted values
    """

    assert pred.shape==truth.shape, "Shapes don't mathc."
    error = np.sqrt(np.mean(np.square(pred - truth), axis=1))
    return error


def get_traj_md2_and_true(location, subjs, abund_type, dict_cv_name, lim,
                          donor="healthy", use_log=True):
    """
    returns the true and mdsine2 predicted trajectory

    (Path) location : directory containing the results of mdsine2 forward sims
    ([str/int]) subjs: List of the subject names
    (str) abund_type : the abundance type
    (dict) dict_cv_name: (str) mdsine2 cv file name -> (str) subject name
    (float) lim : the limit value to be used
    """

    dict_true_abundance = {}
    dict_pred_abundance = {}
    dict_times = {}
    dict_true_abundance_mean = {}
    for sub in subjs:
        cv_name = "{0}-cv{1}".format(donor, sub)
        prefix = dict_cv_name[cv_name]

        true_abund = np.load(location/"{}-cv{}-validate-{}-full-truth"\
                                      ".npy".format(donor, sub, sub))
        pred_abund = np.load(location/"{}-cv{}-validate-{}-full"\
                                      ".npy".format(donor, sub, sub))
        times = np.load(location/"{}-cv{}-validate-{}-full-times"\
                                 ".npy".format(donor, sub, sub))
        if abund_type == "abs":
            true_abund = np.where(true_abund < lim, lim, true_abund)
            true_abund_mean = np.mean(true_abund, axis=1)
            pred_abund = np.where(pred_abund < lim, lim, pred_abund)
            pred_abund_median = np.median(pred_abund, axis=0)
            if use_log:
                dict_true_abundance[prefix] = np.log10(true_abund)
                dict_pred_abundance[prefix] = np.log10(pred_abund_median)
                dict_times[sub] = times
                dict_true_abundance_mean[prefix] = np.log10(true_abund_mean)
            else:
                dict_true_abundance[prefix] = true_abund
                dict_pred_abundance[prefix] = pred_abund_median
                dict_times[sub] = times
                dict_true_abundance_mean[prefix] = true_abund_mean
        elif abund_type == "rel":
            rel_true_abund = true_abund / np.sum(true_abund, axis=0,
                                                 keepdims=True)
            rel_true_abund = np.where(rel_true_abund < lim, lim, rel_true_abund)
            rel_true_abund = rel_true_abund / np.sum(rel_true_abund, axis=0,
                                                     keepdims=True)

            rel_true_abund_mean = np.mean(rel_true_abund, axis=1)
            pred_abund_median = np.median(pred_abund, axis=0)

            rel_pred_abund_median = pred_abund_median / np.sum(pred_abund_median,
                                                               axis=0, keepdims=True)
            rel_pred_abund_median = np.where(rel_pred_abund_median<lim, lim,
                                             rel_pred_abund_median)
            rel_pred_abund_median = rel_pred_abund_median / np.sum(rel_pred_abund_median,
                                                                   axis=0, keepdims=True)

            if use_log:
                dict_true_abundance[prefix] = np.log10(rel_true_abund)
                dict_pred_abundance[prefix] = np.log10(rel_pred_abund_median)
                dict_times[prefix] = times
                dict_true_abundance_mean[prefix] = np.log10(rel_true_abund_mean)
            else:
                dict_true_abundance[prefix] = rel_true_abund
                dict_pred_abundance[prefix] = rel_pred_abund_median
                dict_times[sub] = times
                dict_true_abundance_mean[prefix] = rel_true_abund_mean

    return dict_true_abundance, dict_pred_abundance, dict_true_abundance_mean, dict_times


def get_traj_glv(location, subjs, abund_type, reg, lim, donor="healthy",
                 use_log=True):
    """
    extract the trajectories obtained by running glv model
    (Path) location: location of the file
    ([str]) subjs: list containing the names of the subjects
    (str) abund_type : Relative / Absolute abundance
    (str) reg: type of regression
    (float) the min limit of detection

    @returns
    (dict) (str) subjname name -> (np.ndarray) predicted trajectories
    """

    pred_traj = {}
    for subj in subjs:
        pred_abund = pkl.load(open(location/"glv-{}-{}.pkl".format(reg, subj), "rb")).T
        if abund_type == "abs":
            pred_abund = np.where(pred_abund<lim, lim, pred_abund)
            if use_log:
                pred_traj[subj] = np.log10(pred_abund)
            else:
                pred_traj[subj] = pred_abund
        else:
            rel_pred_abund = pred_abund / np.sum(pred_abund, axis=0,
                keepdims=True)
            rel_pred_abund = np.where(rel_pred_abund<lim, lim, rel_pred_abund)
            rel_true_abund = rel_pred_abund / np.sum(rel_pred_abund, axis=0,
                keepdims=True)
            if use_log:
                pred_traj[subj] = np.log10(rel_pred_abund)
            else:
                pred_traj[subj] = rel_pred_abund

    return pred_traj


def get_traj_clv(location, subjs, model_name, reg_type, lim, use_log=True):

    """
    extracts the trajectories obtained by running relative abundance models
           clv, glv-ra, lra
    (Path) location: location of the file
    ([str]) subjs: list containing the names of the subjects
    (str) model_name : Name of model used
    (str) reg_type: type of regression
    (float) the min limit of detection

    @returns
    (dict) (str) subjname name -> (np.ndarray) predicted trajectories
    """

    pred_traj = {}
    for subj in subjs:
        traj = pkl.load(open(location/"{}-{}-{}.pkl".format(model_name, reg_type,
                                                            subj), "rb")).T
        rel_traj = np.where(traj<lim, lim,  traj)
        rel_traj = rel_traj / np.sum(rel_traj, axis=0, keepdims=True)
        if use_log:
            pred_traj[subj] = np.log10(rel_traj)
        else:
            pred_traj[subj] = rel_traj

    return pred_traj


def combine_into_df(data_dict):
    """
    creates a dataframe contaning data from all subjects and methods
    (dict) data_dict : (str) subject name -> (dict) (str) method -> (np.ndarray)
    """

    all_data_y = {}
    for key1 in data_dict:
        for key2 in data_dict[key1]:
            if key2 not in all_data_y:
                all_data_y[key2] = []
            all_data_y[key2] += list(data_dict[key1][key2])

    methods = []
    errors = []
    for keys in all_data_y:
        n = len(all_data_y[keys])
        methods += [keys] * n
        errors += all_data_y[keys]

    df = pd.DataFrame(list(zip(methods, errors)), columns=["Method", "Error"])
    return df

def format_box_plot(pred_dict, true_dict, type_):
    """
    (dict) pred_dict: (str) method_name -> (dict) (str) subject name ->(np.ndarray)
    (dict) true_dict: (str) subj name -> (np.ndarray) the true values
    """
    error_dict = {}

    for k1 in pred_dict:
        for k2 in pred_dict[k1]:
            if k2 not in error_dict:
                error_dict[k2] = {}
            error_dict[k2][k1] = compute_rms_error(pred_dict[k1][k2], true_dict[k2])
    combined_df = combine_into_df(error_dict)
    test_results = signed_rank_test_box(error_dict, "MDSINE2", type_)

    return combined_df, test_results

def signed_rank_test_box(data_dict, ref_key, type_):
    """
    performs Wilcoxon signed rank test
    (dict) data_dict
    (str) ref_key: one of the keys (method_names) in data_dict with which the
    comparison is done
    """
    # print("running Wilcoxon signed-rank test for {}".format(donor))
    combined_data = {}
    for key1 in data_dict:
        for key2 in data_dict[key1]:
            if key2 not in combined_data:
                combined_data[key2] = []
                # print(key1, key2)
            combined_data[key2] += list(data_dict[key1][key2])

    ref_data = combined_data[ref_key]
    p_vals = []
    to_use_order = []
    order = []
    if type_ =="rel":
        to_use_order = REL_ORDER
    elif type_ =="abs":
        to_use_order = ABS_ORDER

    for key in to_use_order:
        if key != ref_key:
            order.append("{} and {}".format(ref_key, key))
            s, p = stats.wilcoxon(ref_data, combined_data[key], alternative='less')
            p_vals.append(p)

    test = multitest.multipletests(p_vals, alpha=0.05, method="fdr_bh")

    return list(test[1])


def combine_data_alternative(pred_dict, mean_dict, binned_dict):
    """combines the errors and mean abundances in a data frame used
       for the alternative plot"""
    pred_method_dict = {}
    mean_method_dict = {}
    bin_method_dict = {}

    for key1 in pred_dict:
        for key2 in pred_dict[key1]:
            if key2 not in pred_method_dict:
                pred_method_dict[key2] = []
            if key2 not in mean_method_dict:
                mean_method_dict[key2] = []
            if key2 not in bin_method_dict:
                bin_method_dict[key2] = []
            pred_method_dict[key2] += list(pred_dict[key1][key2])
            mean_method_dict[key2] += list(mean_dict[key1])
            bin_method_dict[key2] += list(binned_dict[key1])

    bins = []
    methods = []
    errors = []
    mean_abund = []

    for keys in pred_method_dict:
        n = len(pred_method_dict[key2])
        methods += [keys] * n
        errors += pred_method_dict[keys]
        mean_abund += mean_method_dict[keys]
        bins += bin_method_dict[keys]

    df = pd.DataFrame(list(zip(methods, errors, bins, mean_abund)),
                      columns=["Method", "Error", "Floor", "Mean Abundance"])

    return df


def format_alternative_plot(pred_dict, true_dict, mean_abund_dict):
    """
    formats the data in a format suitable for the alternative plot which breaks down
    errors for different abundance level

    (dict) pred_dict: (str) method_name -> (dict)(str) subject name ->
                       (np.ndarray) 2-dim array containing predicted abundances
    (dict) true_dict: (str) subject name -> (np.ndarray) 2-dim array containing
                       absolute abundances
    (dict) mean_abund_dict: (str) subject name -> (np.ndarray) mean abundances
    """

    error_dict = {}
    abund_all = []
    binned_data = {}
    for k1 in pred_dict:
        for k2 in sorted(list(pred_dict[k1].keys())):
            if k2 not in error_dict:
                error_dict[k2] = {}
            error_dict[k2][k1] = compute_rms_error(pred_dict[k1][k2], true_dict[k2])

    for k in sorted(list(mean_abund_dict.keys())):
        abund_all += list(mean_abund_dict[k])
        binned_data[k] = np.floor(mean_abund_dict[k])

    bar_df = pd.DataFrame(abund_all, columns=["Mean Abundance"])
    combined_df = combine_data_alternative(error_dict, mean_abund_dict, binned_data)

    return combined_df, bar_df


def add_annotation(axes, data_df, x_var, y_var, box_pairs, p_values, order):
    """add annotations to axis """
    def star_p(p):
        if p >=0.05:
            return "ns"
        elif 0.01 <= p < 0.05:
            return "*"
        elif 0.001 <= p < 0.01:
            return "**"
        elif 0.0001 <= p < 0.001:
            return "***"
        else:
            return "****"

    order_index = {order[i]:i for i in range(len(order))}
    y_prev = 0
    for i in range(len(box_pairs)):
        pair = box_pairs[i]
        y1_data = data_df.loc[data_df[x_var]==pair[0]]
        y1 = np.percentile(y1_data[y_var].to_numpy(), 97.5)
        y2_data = data_df.loc[data_df[x_var]==pair[1]]
        y2 = np.percentile(y2_data[y_var].to_numpy(), 97.5)
        y = max(y1, y2)+0.025
        if y_prev ==0:
            y_prev = y
        if y < y_prev:
            y = y_prev
        x1 = order_index[pair[0]]
        x2 = order_index[pair[1]]
        h = 0.1
        axes.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c="k")
        star = star_p(p_values[i])
        axes.text((x1+x2)/2, y + h, star, ha = "center", va="bottom", color="k")
        y_prev = y+ 3*h


def box_plot(data_df, axes, title, use_log, type_, ylab, pvalues):
    """
    make a box plot
    (pd.DataFrame) data_df : df containing the data to be plotted
    (bool) use_log: whether or not to use log scale
    ([float]) pvalues: list of p-values obtained from signed rank-test
    """

    axes.xaxis.grid(False)
    axes.yaxis.grid(True)

    axes.set_title(title, loc="left", fontweight="bold", fontsize=TITLE_FONTSIZE)
    if not use_log:
        axes.set_yscale("log")
    palette = PAL_REL
    order = REL_ORDER
    if type_ == "abs":
        palette = PAL_ABS
        order = ABS_ORDER

    sns.boxplot(y="Error", x="Method", data=data_df, whis=[2.5, 97.5], width=.75,
                showfliers=False, ax=axes, palette=palette, order=order)
    sns.stripplot(y="Error", x="Method", data=data_df, size=2,
                  linewidth=0.5, alpha=0.5, ax=axes, palette=palette, order=order) #, color=".3"

    axes.set_ylabel(ylab, fontsize=AXES_FONTSIZE, fontweight="bold")
    axes.set_xlabel("Model", fontsize=AXES_FONTSIZE, labelpad=3, fontweight="bold")
    axes.set_xticklabels(axes.get_xticklabels(), rotation=0, fontsize=TICK_FONTSIZE)
    axes.tick_params(axis='y', labelsize=TICK_FONTSIZE)

    if type_=="abs":
        add_annotation(axes, data_df, "Method", "Error",
        [("MDSINE2", "gLV-ridge"), ("MDSINE2", "gLV-elastic\n net")], pvalues,
                       ABS_ORDER)
        axes.set_ylim([0, 4])
    else:
        add_annotation(axes, data_df, "Method", "Error",
                       [("MDSINE2","cLV"),("MDSINE2","LRA"), ("MDSINE2","gLV-RA"),
                       ("MDSINE2", "gLV-ridge"),
                       ("MDSINE2", "gLV-elastic\n net")], pvalues,
                       REL_ORDER)
        axes.set_ylim([0, 5])


def get_bin_category(bins, df):

    df_np = df["Mean Abundance"].to_numpy()
    bins_np = np.asarray(bins)

    bins_ = np.digitize(df_np, bins_np)
    final_bins = []
    for i in bins_:
        #if i==5:
            #print("yes")
        if i == len(bins):
            final_bins.append(len(bins)-1)
        else:
            final_bins.append(i)

    return np.asarray(final_bins)


def alternative_plot(data_df, title1, title2, use_log, bar_df, type_, axes1,
                     axes2, y_lab, donor, add_legend=False):
    """makes the alternative plot based on the DataFrame data_df"""

    def test_stars(p, dtype):
        star = ""
        if p < 0.05:
            star= "x"#"$\\times$"
        else:
            star= "o"#"$\\circ$"
        #print(p, star)
        return star

    quantile_10 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    bins = bar_df["Mean Abundance"].quantile(quantile_10).to_list()
    bin_category = get_bin_category(bins, data_df)
    data_df["Bins"] = bin_category
    order, raw_p, p_values = signed_rank_test_alt(data_df, donor, type_)


    xtick_labels = ["[{:.1f},\n {:.1f}]".format(bins[i-1], bins[i])
                    for i in range(1, len(bins))]

    if type_ =="rel":
        xtick_labels = ["[{:.1f},\n {:.1f}]".format(bins[i-1], bins[i])
                        for i in range(1, len(bins))]

    #export_p(order, raw_p, p_values, xtick_labels, donor+"_"+type_)
    palette = PAL_REL
    order = REL_ORDER
    if type_ =="abs":
        palette = PAL_ABS
        order = ABS_ORDER
    sns.boxplot(y="Error", x="Bins", hue="Method", data=data_df, whis=[2.5, 97.5],
                showfliers=False, ax=axes2, palette=palette, hue_order=order)#, palette=my_pal
    sns.stripplot(y="Error", x="Bins", hue="Method", data=data_df, size=2,
                  alpha=0.5, ax=axes2, dodge=True, palette=palette,
                  linewidth=0.5, hue_order=order) #color=".8",

    y_data = data_df.groupby(["Bins", "Method"])["Error"].max().values
    y_percentile = data_df.groupby(["Bins", "Method"])["Error"].quantile(1).values
    star = [test_stars(p, type_) for p in p_values]
    y_data_round = ["{:.1f}".format(p) for p in y_data]

    handles, labels = axes2.get_legend_handles_labels()
    histplot = sns.histplot(data=bar_df.to_numpy(), ax=axes1, bins=bins,
                            color="lime", cbar=False, legend=False, stat="density")

    i_star = 0
    i_box = 0

    for tick in range(len(axes2.get_xticklabels())):
        if type_ == "abs":
            y = np.max(y_percentile[i_box+1:i_box+3]) +0.05
            if y > 4.5:
                y = 4.6
            axes2.text(tick+0.145, y, "   ".join(star[i_star:i_star+2]),
                       horizontalalignment="center", color="black", fontsize=10)
        else:
            y = np.max(y_percentile[i_box+1:i_box+6])
            if y > 4.5:
                y = 4.5

            axes2.text(tick+0.045, y, "  ".join(star[i_star:i_star+5]),
                       horizontalalignment="center", color="black", fontsize=10)

        if type_ =="abs":
            i_star += 2
            i_box += 3
        elif type_=="rel":
            i_star += 5
            i_box += 6

    if type_ =="rel":
        axes1.set_xlabel("log (Mean Rel Abundance)", fontsize=AXES_FONTSIZE,
                         labelpad=3, fontweight="bold")
        axes2.set_xlabel("log (Mean Rel Abundance)", fontsize=AXES_FONTSIZE,
                         fontweight="bold")
        axes1.set_ylim(0, 0.5)
        axes2.set_ylim([0, 5])

    else:
        axes1.set_xlabel("log (Mean Abs Abundance)", fontsize=AXES_FONTSIZE,
                         labelpad=3, fontweight="bold")
        axes2.set_xlabel("log (Mean Abs Abundance)", fontsize=AXES_FONTSIZE,
                         fontweight="bold")
        axes1.set_ylim(0, 0.5)
        axes2.set_ylim([0, 4])

    axes1.set_title(title1, loc="left", fontweight="bold", fontsize=TITLE_FONTSIZE)
    axes1.set_ylabel("Frequency", fontsize=AXES_FONTSIZE,
                     fontweight="bold")
    axes2.set_ylabel(y_lab, fontsize=AXES_FONTSIZE, fontweight="bold")

    axes1.tick_params(axis='y', labelsize=TICK_FONTSIZE)
    axes1.tick_params(axis='x', labelsize=TICK_FONTSIZE)
    axes2.tick_params(axis='y', labelsize=TICK_FONTSIZE)

    axes2.yaxis.grid(True)
    axes2.set_title(title2, loc="left", fontweight="bold", fontsize=TITLE_FONTSIZE)
    axes2.set_xticklabels(xtick_labels, fontsize=TICK_FONTSIZE)

    if not use_log:
        axes2.set_yscale("log")

    if add_legend:
        n = 3
        if type_ =="rel":
            n = 6
        lgd2 = axes2.legend(handles[0:n], labels[0:n], bbox_to_anchor=(1.01, 1),
                          loc=2, borderaxespad=0., fontsize=TICK_FONTSIZE,
                          title_fontsize=AXES_FONTSIZE, title="$\\bf{Model}$")
        axes2.add_artist(lgd2)

        handles = []
        l = mlines.Line2D([],[], color="black", linestyle='none',
                          marker='x', label='p < 0.05', markerfacecolor='none')
        handles.append(l)
        l = mlines.Line2D([],[], color="black",linestyle='none',
            marker ='o', label='p > 0.05', markerfacecolor='none')
        handles.append(l)
        lgnd3 = axes2.legend(handles = handles, title='$\\bf{P-values}$',
                             loc='lower left', borderaxespad=0., bbox_to_anchor=(1.01, 0),
                             title_fontsize=TICK_FONTSIZE, fontsize = TICK_FONTSIZE)
        axes2.add_artist(lgnd3)
    else:
        axes2.get_legend().remove()

def signed_rank_test_alt(data_df, donor, type_):

    bins_dict = {}
    otu_order = {}
    for row in data_df.to_numpy():
        method = row[0]
        bin_id = row[-1]
        mean_abundance = row[3]
        error = row[1]
        if bin_id not in bins_dict:
            bins_dict[bin_id] = {}
        if bin_id not in otu_order:
            otu_order[bin_id] = {}
        if method not in bins_dict[bin_id]:
            bins_dict[bin_id][method] = []
        if method not in otu_order[bin_id]:
            otu_order[bin_id][method] = []

        bins_dict[bin_id][method].append(error)
        otu_order[bin_id][method].append(mean_abundance)

    ref_key = "MDSINE2"
    p_values = []
    order  = []
    to_use_order = ""
    if type_ =="rel":
        to_use_order = REL_ORDER
    elif type_ =="abs":
        to_use_order = ABS_ORDER

    for b in sorted(bins_dict.keys()):
        ref_data = bins_dict[b][ref_key]
        #print("bin:", b)
        for m in to_use_order:
            if m != ref_key:
                other_data = bins_dict[b][m]
                order.append("{} and {}".format(ref_key, m))
                s, p = stats.wilcoxon(ref_data, other_data, alternative='less')
                p_values.append(p)
    test = multitest.multipletests(p_values, alpha=0.05, method="fdr_bh")

    return order, p_values, test[1]

def plot(true_val, pred_vals, x_data, taxa_info, dir, mouse_num, mouse_study,
         y_lab="Abundance (CFUs/g)"):

    N = true_val.shape[0]
    dir.mkdir(parents=True, exist_ok=True)

    for i in range(N):
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)
        axes.set_xlabel("Times")
        axes.set_ylabel(y_lab)
        for key in pred_vals:
            axes.plot(x_data, 10**pred_vals[key][i], label=key)
        axes.plot(x_data, 10**true_val[i], label="true")
        axes.legend()
        axes.set_yscale("log")
        name = "{}, {}, mouse {}".format(taxa_info[i].name, md2.pylab.base.taxaname_for_paper(
            taxa_info[i], taxa_info), mouse_num)
        axes.set_title(name)
        md2.visualization.shade_in_perturbations(axes, perturbations=mouse_study.perturbations,
                                                 subj=mouse_study, textsize=8)

        fig.savefig(dir/"{}.pdf".format(taxa_info[i].name))
        plt.close()


def trajectory_plotter(true_dict, pred_dict, destination, study_file_loc, subjects,
                       dict_subj, ylabel, use_log=True):
    """
    plots the trajectories for different methods
    (dict) true_dict: (str) subj_name -> (np.ndarray)
    (dict) pred_dict:(str) method_name -> (dict) (str) subj_name -> (np.ndarray)
    (Path) destination: directory where the plots are saved
    (Path) study_file_loc: pickle file containing the location
    """

    pred_dict_key_swapped = {}
    for key1 in pred_dict:
        for key2 in pred_dict[key1]:
            if key2 not in pred_dict_key_swapped:
                pred_dict_key_swapped[key2] = {}
            pred_dict_key_swapped[key2][key1] = pred_dict[key1][key2]


    study_file = md2.Study.load(study_file_loc)
    taxa = study_file.taxa
    destination = Path(destination)
    for sub in subjects:
        true_traj = true_dict[sub]
        subj_pkl = study_file[dict_subj[sub]]
        times = subj_pkl.times
        plot(true_traj, pred_dict_key_swapped[sub], times, taxa,
             destination/"mouse_{}".format(dict_subj[sub]), dict_subj[sub],
             study_file[dict_subj[sub]], y_lab=ylabel)


def main():
    print("Making Figure 3")
    args = parse_args()
    healthy_subjs = ["2", "3", "4", "5"]
    dict_subj = {"0":"2", "1":"3", "2":"4", "3":"5"}
    use_log = True
    rel_lim = args.rel_abund_lim
    abs_lim = 1e5

    fig = plt.figure(figsize=(24, 13))
    spec = gridspec.GridSpec(ncols=40, nrows=21, figure=fig)

    #for the box plots showing the error
    ax_he_abs_box = fig.add_subplot(spec[0:6, 0:14])
    ax_he_rel_box = fig.add_subplot(spec[0:6, 16:40])
    #for the histograms and box plots
    ax_he_abs_hist = fig.add_subplot(spec[8:11, 0:14])
    ax_he_abs_bin = fig.add_subplot(spec[13:21, 0:14])
    ax_he_rel_hist = fig.add_subplot(spec[8:11, 16:40])
    ax_he_rel_bin = fig.add_subplot(spec[13:21, 16:40])

    dict_cv = {"healthy-cv5": "3", "healthy-cv4": "2", "healthy-cv3": "1",
               "healthy-cv2": "0"}

    healthy_subjs_clv = [dict_cv[key] for key in dict_cv]

    mdsine_path = args.mdsine_path
    clv_elas_path = args.non_glv_elas_path
    clv_ridge_path = args.non_glv_ridge_path
    glv_elas_path = args.glv_elas_path
    glv_ridge_path = args.glv_ridge_path

    true_abund_abs, pred_abund_md2_median_abs, true_abund_mean_abs, times = get_traj_md2_and_true(
        Path(args.mdsine_path), healthy_subjs, "abs", dict_cv, abs_lim)
    true_abund_rel, pred_abund_md2_median_rel, true_abund_mean_rel, times = get_traj_md2_and_true(
        Path(args.mdsine_path), healthy_subjs, "rel", dict_cv, rel_lim)

    pred_glv_ridge_abs = get_traj_glv(Path(glv_ridge_path), healthy_subjs_clv,
        "abs", "ridge", abs_lim)
    pred_glv_elastic_abs = get_traj_glv(Path(glv_elas_path), healthy_subjs_clv,
        "abs", "elastic-net", abs_lim)

    pred_glv_ridge_rel = get_traj_glv(Path(glv_ridge_path), healthy_subjs_clv,
        "rel", "ridge", rel_lim)
    pred_glv_elastic_rel = get_traj_glv(Path(glv_elas_path), healthy_subjs_clv,
        "rel", "elastic-net", rel_lim)

    pred_clv_rel = get_traj_clv(Path(clv_elas_path), healthy_subjs_clv, "clv",
                                "elastic-net", rel_lim)
    pred_lra_rel = get_traj_clv(Path(clv_elas_path), healthy_subjs_clv, "lra",
                                "elastic-net", rel_lim)
    pred_glv_ra = get_traj_clv(Path(clv_elas_path), healthy_subjs_clv, "glv-ra",
                                "elastic-net", rel_lim)
    #pred_glv_ridge_rel = get_traj_clv(Path(clv_ridge_path), healthy_subjs_clv, "glv",
    #                            "ridge", rel_lim)
    #pred_glv_elastic_rel = get_traj_clv(Path(clv_elas_path), healthy_subjs_clv, "glv",
    #                            "elastic-net", rel_lim)


    dict_pred_abs = {"MDSINE2":pred_abund_md2_median_abs,"gLV-ridge": pred_glv_ridge_abs,
                     "gLV-elastic\n net": pred_glv_elastic_abs}
    dict_pred_rel = {"MDSINE2":pred_abund_md2_median_rel, "cLV": pred_clv_rel,
                     "gLV-elastic\n net": pred_glv_elastic_rel, "LRA": pred_lra_rel,
                     "gLV-RA": pred_glv_ra, "gLV-ridge": pred_glv_ridge_rel}

    healthy_abs_box_df, test_healthy_abs = format_box_plot(dict_pred_abs,
                                                          true_abund_abs, "abs")
    healthy_rel_box_df, test_healthy_rel = format_box_plot(dict_pred_rel,
                                                           true_abund_rel, "rel")
    healthy_rel_box_df.to_csv("healthy_rel.csv", sep=",")

    healthy_abs_alt_df, healthy_abs_abund_df = format_alternative_plot(dict_pred_abs,
                                                                       true_abund_abs,
                                                                       true_abund_mean_abs)
    healthy_rel_alt_df, healthy_rel_abund_df = format_alternative_plot(dict_pred_rel,
                                                                       true_abund_rel,
                                                                       true_abund_mean_rel)
    fig = plt.figure(figsize=(24, 13))
    spec = gridspec.GridSpec(ncols=40, nrows=21, figure=fig)

    ax_he_abs_box = fig.add_subplot(spec[0:6, 0:14])
    ax_he_rel_box = fig.add_subplot(spec[0:6, 16:40])

    ax_he_abs_hist = fig.add_subplot(spec[8:11, 0:14])
    ax_he_abs_bin = fig.add_subplot(spec[13:21, 0:14])

    ax_he_rel_hist = fig.add_subplot(spec[8:11, 16:40])
    ax_he_rel_bin = fig.add_subplot(spec[13:21, 16:40])

    box_plot(healthy_abs_box_df, ax_he_abs_box, "A", use_log, "abs",
        "RMSE (log Abs Abundance)", test_healthy_abs)
    box_plot(healthy_rel_box_df, ax_he_rel_box, "B", use_log, "rel",
        "RMSE (log Rel Abundance)", test_healthy_rel)

    alternative_plot(healthy_abs_alt_df, "C", "E", use_log, healthy_abs_abund_df,
       "abs", ax_he_abs_hist, ax_he_abs_bin, "RMSE (log Abs Abundance)",
       "healthy", add_legend=False)

    alternative_plot(healthy_rel_alt_df, "D", "F", use_log, healthy_rel_abund_df,
       "rel", ax_he_rel_hist, ax_he_rel_bin, "RMSE (log Rel Abundance)",
        "healthy", add_legend=True)

    output_path = Path(args.output_path)
    output_path.mkdir(exist_ok=True, parents=True)
    fig.savefig(output_path/"{}.pdf".format(args.output_name), bbox_inches="tight")

    # uncomment to plot the trajectories
    #trajectory_plotter(true_abund_abs, dict_pred_abs,
    #    "trajectories_absolute", args.study_file,
    #    healthy_subjs_clv, dict_subj, "CFUs/g")

    #trajectory_plotter(true_abund_rel, dict_pred_rel,
    #    "trajectories_relative", args.study_file,
    #    healthy_subjs_clv, dict_subj, "Relative abundance")

main()
