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
#print(font_files)

for font_file in font_files:
    #print(font_file)
    ff = font_file.split("/")[-1]
    if "._" not in ff:
        font_manager.fontManager.addfont(font_file)

# change font
rcParams['font.family'] = 'Arial'

def parse_args():

    parser = argparse.ArgumentParser(description = "files needed for making"\
    "main figure 2")
    parser.add_argument("-loc", "--deseq_loc", required = "True",
        help = "location of the folder containing the results of deseq")
    parser.add_argument("-abund", "--abundance", required='True',
        help="High/Low abundance")
    parser.add_argument("-txt", "--txt_file", required = "True",
        help = "name of the text file containing the names of orders")
    parser.add_argument("-taxo", "--taxonomy", required = "True",
        help = "name of the taxonomy hierarchy")
    parser.add_argument("-o", "--output_name", required = "True",
        help = "Name of the output file")
    parser.add_argument("-o_loc", "--output_loc", required="True",
        help = "location of the output")

    return parser.parse_args()


def get_names(df):

    names = []
    index = df.index
    np_df = df.to_numpy()
    for i in range(np_df.shape[0]):
        p = np_df[i, -1]
        if p <= 0.05:
            names.append(index[i])

    return set(names)


def get_deseq_info(loc, donor):
    """
    extracts relevant information from the csv file containing the deseq results
    """

    hfd_names = get_names(pd.read_csv(loc / "{}1.csv".format(donor),
                         index_col=0))
    vanc_names = get_names(pd.read_csv(loc /"{}2.csv".format(donor),
                         index_col=0))
    gent_names = get_names(pd.read_csv(loc / "{}3.csv".format(donor),
                         index_col=0))

    names_union = hfd_names.union(vanc_names.union(gent_names))
    cleaned_names = []
    names_dict = {}
    i = 0

    for name in names_union:
        name_split = name.split()
        order_family = " ".join(name_split[2:])
        class_ = name_split[1]
        phylum = name_split[0]
        if "unknown" in order_family:
            if order_family !="unknown unknown":
                cleaned_names.append(order_family.replace("unknown", "unknown"))
            else:
                if class_ == "unknown":
                    cleaned_names.append("{},phylum".format(phylum) + " " +
                            "Phylum")
                else:
                    cleaned_names.append("{},class".format(class_) + " " +
                            "Class")
        else:
            cleaned_names.append(order_family)
        names_dict[cleaned_names[-1]] = name
        i += 1

    return set(cleaned_names), names_dict

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

def make_df(loc, names, dict_, donor):
    """generates the df used to make the heatmap"""

    hfd_df = pd.read_csv(loc / "{}1.csv".format(donor),
                         index_col=0)
    vanc_df = pd.read_csv(loc / "{}2.csv".format(donor),
                         index_col=0)
    gent_df = pd.read_csv(loc / "{}3.csv".format(donor),
                         index_col=0)

    dict_df = {"High Fat Diet": [], "Vancomycin":[], "Gentamicin":[]}

    dict_df["High Fat Diet"] = select_rows(names, hfd_df, dict_)
    dict_df["Vancomycin"] = select_rows(names, vanc_df, dict_)
    dict_df["Gentamicin"] = select_rows(names, gent_df, dict_)

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

    fig = plt.figure(figsize=(3.5, 11))
    gs = fig.add_gridspec(28, 6)
    axes1 = ""
    axes2 = ""

    if taxonomy == "order":
        if abundance == "high":
            axes1 = fig.add_subplot(gs[0:18,  0:3])
        elif abundance == "low":
            axes1 = fig.add_subplot(gs[19:28, 0:3])

    if taxonomy == "phylum":
        gs = fig.add_gridspec(28, 7)
        if abundance == "high":
            axes1 = fig.add_subplot(gs[0:18,  0:3])
        elif abundance == "low":
            axes1 = fig.add_subplot(gs[19:28, 0:3])

    make_heatmap(df1, axes1, False, "Healthy", False)

    loc.mkdir(exist_ok=True, parents=True)
    fig.savefig(loc / "{}.pdf".format(name), bbox_inches="tight")
    plt.close()

def change_labels(df1, df2):
    """adjust the labels in (pd.DataFrame) df1 and df2"""

    def find_lengths(index1, index2):
        max_len_order = 0
        max_string_len = 0
        combined = list(index1) + list(index2)
        refined_names_dict = {}
        for name in combined:
            refined_name = name
            if "_" in name:
                refined_name = "_".join(name.split("_")[0:-1])
            if "unknown" in name and "phylum" not in name:
                #print("unknown:", name)
                refined_name = name.replace("unknown", "(not resolved to family)")
            if "phylum" in name:
                refined_name = "(" + name.split(",")[0] + ")" + " (not resolved"+\
                    " below phylum)"
            len_order = len(refined_name.split()[0])
            refined_names_dict[name] = refined_name.replace("_", " ")
            if len_order > max_len_order:
                max_len_order = len_order
            if len(refined_name) > max_string_len:
                max_string_len = len(refined_name)
        return max_len_order, max_string_len -max_len_order, refined_names_dict

    def get_final(N_order, N_total, names, names_dict):
        final_names = []
        for name in names:
            split_name = names_dict[name].split()
            d1 = N_order -len(split_name[0]) + 2
            d2 = N_total - len(" ".join(split_name[1:])) + 2
            new_name = split_name[0] + " " * 1 + " ".join(split_name[1:]) + " "*1
            final_names.append(new_name)

        return final_names

    index1 = df1.index
    index2 = df2.index
    order_len, string_len, dict_ = find_lengths(index1, index2)
    final_name1 = get_final(order_len, string_len, index1, dict_ )
    final_name2 = get_final(order_len, string_len, index2, dict_ )
    df1.index = final_name1
    df2.index = final_name2

    return df1, df2


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
            linecolor="black", linewidths=1, ax=axes, vmin=-2, vmax=2,
            xticklabels=label_x)
        axes.set_xticklabels(axes.get_xticklabels(), rotation=75, fontsize=12)
    else:
        axes = sns.heatmap(df, cmap=cmap, cbar=False, yticklabels=label_y,
            linecolor="black", linewidths=1, ax=axes, xticklabels=False, vmin=-2,
            vmax=2)

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

if __name__ =="__main__":

    print("Making Heatmap Order")
    args = parse_args()
    deseq_loc = Path(args.deseq_loc)
    file = open(deseq_loc /"{}.txt".format(args.txt_file))
    labels = file.read().split(", ")
    new_labels = []
    for lab in labels:
        new_lab = lab
        if args.taxonomy == "phylum":
            new_lab = lab.split()[1]
        if "NA" in lab:
            new_lab = lab.replace("NA", "unknown")
        new_labels.append(new_lab)

    healthy_names, healthy_dict = get_deseq_info(deseq_loc, "healthy_order")
    names_union = healthy_names.union(healthy_names)

    all_dict = {}
    for keys in healthy_dict:
        if keys not in all_dict:
            all_dict[keys] = healthy_dict[keys]
    names_not_abundant = sorted(list(set(names_union) - set(new_labels)))
    healthy_df_abundant = make_df(deseq_loc, new_labels, all_dict, "healthy_order")

    healthy_df_non_abundant = make_df(deseq_loc, names_not_abundant, all_dict,
        "healthy_order")

    healthy_df_abundant, healthy_df_non_abundant = change_labels(
        healthy_df_abundant, healthy_df_non_abundant)


    if args.abundance == "high":
        make_plot(healthy_df_abundant, args.abundance,
            args.taxonomy, args.output_name, Path(args.output_loc))
    elif args.abundance == "low":
        make_plot(healthy_df_non_abundant, args.abundance,
            args.taxonomy, args.output_name, Path(args.output_loc))
