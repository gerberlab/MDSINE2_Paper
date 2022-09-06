import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os

from matplotlib import rcParams
from matplotlib import font_manager
from pathlib import Path

rcParams['pdf.fonttype'] = 42

font_dirs = ['figures/arial_fonts']
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

for font_file in font_files:
    ff = font_file.split("/")[-1]
    if "._" not in ff:
        font_manager.fontManager.addfont(font_file)

# change font
rcParams['font.family'] = 'Arial'


def select_rows(names, df, dict_):
    """selects rows associated with name in names from the df (DataFrame)"""

    values = []
    for short_name in names:
        if short_name in dict_:
            name = dict_[short_name]
            if name in df.index:
                row = df.loc[name]
                p = row["padj"]
                log_change = row["log2FoldChange"]
                if p <=0.05:
                    if log_change < 0:
                        values.append(-1)
                    else:
                        values.append(1)
                else:
                    values.append(np.nan)
            else:
                values.append(np.nan)
        else:
            values.append(np.nan)
    return values


def make_df(files_loc, file_name, names, names_dict):
    """genetates the dataframes used to make the heatmap"""

    hfd_df = pd.read_csv(files_loc / "{}1.csv".format(file_name),
                         index_col=0)
    vanc_df = pd.read_csv(files_loc / "{}2.csv".format(file_name),
                         index_col=0)
    gent_df = pd.read_csv(files_loc / "{}3.csv".format(file_name),
                         index_col=0)
    dict_df = {"High Fat Diet": [], "Vancomycin":[], "Gentamicin":[]}

    dict_df["High Fat Diet"] = select_rows(names, hfd_df, names_dict)
    dict_df["Vancomycin"] = select_rows(names, vanc_df, names_dict)
    dict_df["Gentamicin"] = select_rows(names, gent_df, names_dict)

    combined_df = pd.DataFrame.from_dict(dict_df)
    combined_df.index = names

    return combined_df


def make_plot(df1, abundance, taxonomy, name, loc):
    """
    sets up the parameters needed to make the heatmap
    (pd.DataFrame) df1: the data frame to be plotted
    (str) abundance : high / low abundance
    (str) taxonomy: the taxonomy hierarchy under consideration
    (str) name : the name to save the plot
    (Path) loc: directory where the figure is saved
    """

    def make_heatmap(df, axes, label_x=False, title=None, label_y=False):
        """
        makes the heatmap showing deseq results
        (pd.DataFrame) df: the dataframe containing the deseq results
        (plt.axis) axes : the axes on which the heatmap is plotted
        (bool) label_x, label_y : whether or not to label the axes
        (str) title: title of the plot
        """

        colors = {"red":-1, "blue":1}
        l_colors = sorted(colors, key=colors.get)
        #cmap = mpl.colors.ListedColormap(l_colors)
        cmap = "RdBu"
        axes.set_aspect("equal")
        if label_x:
            axes = sns.heatmap(df, cmap=cmap, cbar=False, yticklabels=label_y,
                               linecolor="black", linewidths=1, ax=axes, vmin=-2,
                               vmax=2, xticklabels=label_x)
            axes.set_xticklabels(axes.get_xticklabels(), rotation=75, fontsize=12)
        else:
            axes = sns.heatmap(df, cmap=cmap, cbar=False, yticklabels=label_y,
                               linecolor="black", linewidths=1, ax=axes,
                               xticklabels=False, vmin=-2, vmax=2)

    #if title is not None:
    #    axes.set_title(title, fontweight="bold", fontsize=16)

        data = df.to_numpy()
        for y in range(data.shape[0]):
            for x in range(data.shape[1]):
                n = data[y, x]
                text = ""
                if not np.isnan(n):
                    if n > 0:
                        axes.text(x + 0.5, y+0.5, "$+$", horizontalalignment='center',
                                  verticalalignment='center', fontsize=12)
                    else:
                        axes.text(x + 0.5, y+0.5, "$-$", horizontalalignment='center',
                                  verticalalignment='center', fontsize=12)
        for _, spine in axes.spines.items():
            spine.set_visible(True)
    fig = plt.figure(figsize=(3.5, 11))
    gs = fig.add_gridspec(28, 6)
    axes1 = ""
    axes2 = ""

    if taxonomy == "order":
        if abundance == "high_abund_defined":
            axes1 = fig.add_subplot(gs[0:18,  0:3])
        elif abundance == "high_abund_not_defined":
            axes1 = fig.add_subplot(gs[19:23,  0:3])
        elif abundance == "low_abund":
            axes1 = fig.add_subplot(gs[23:28, 0:3])


    if taxonomy == "phylum":
        gs = fig.add_gridspec(28, 7)
        if abundance == "high":
            axes1 = fig.add_subplot(gs[0:18,  0:3])
        elif abundance == "low":
            axes1 = fig.add_subplot(gs[19:28, 0:3])

    make_heatmap(df1, axes1, False, "Healthy", True)

    loc.mkdir(exist_ok=True, parents=True)
    fig.savefig(loc / "{}.pdf".format(name), bbox_inches="tight")
    plt.close()
