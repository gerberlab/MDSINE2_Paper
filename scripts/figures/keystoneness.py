#@title
import numpy as np
import mdsine2 as md2
from tqdm.notebook import tqdm
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm, gridspec
import matplotlib.colors as mcolors
import seaborn as sns
from mdsine2.names import STRNAMES
import scipy
import scipy.stats


def cluster_nonmembership_df(md):
    entries = []
    for cluster in md.get_clustering():
        for otu in md.taxa:
            if otu.idx not in cluster.members:
                entries.append({
                    "ClusterID": cluster.id,
                    "OTU": otu.name
                })
    return pd.DataFrame(entries)


def cluster_membership_df(md):
    entries = []
    for cluster in md.get_clustering():
        for oidx in cluster.members:
            otu = md.get_taxa(oidx)
            entries.append({
                "ClusterOfOTU": cluster.id,
                "OTU": otu.name
            })
    return pd.DataFrame(entries)


def create_cmap(tag, nan_value="red"):
    cmap = cm.get_cmap(tag)
    cmap.set_bad(color=nan_value)
    return cmap


# Preprocessing for Cycles (Figure F)
class MdsineOutput(object):
    '''
    A class to encode the data output by MDSINE.
    '''
    def __init__(self, dataset_name, pkl_path):
        self.dataset_name = dataset_name
        self.mcmc = md2.BaseMCMC.load(pkl_path)
        self.taxa = self.mcmc.graph.data.taxa
        self.name_to_taxa = {otu.name: otu for otu in self.taxa}

        self.interactions = None
        self.clustering = None

        self.clusters_by_idx = {
            (c_idx): [self.get_taxa(oidx) for oidx in cluster.members]
            for c_idx, cluster in enumerate(self.get_clustering())
        }

    @property
    def num_samples(self) -> int:
        return self.mcmc.n_samples

    def get_cluster_df(self):
        return pd.DataFrame([
            {
                "id": cluster.id,
                "idx": c_idx+1,
                "otus": ",".join([self.get_taxa(otu_idx).name for otu_idx in cluster.members]),
                "size": len(cluster)
            }
            for c_idx, cluster in enumerate(self.clustering)
        ])

    def get_interactions(self):
        if self.interactions is None:
            self.interactions = self.mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk(section='posterior')
        return self.interactions

    def get_taxa(self, idx):
        return self.taxa.index[idx]

    def get_taxa_by_name(self, name: str):
        return self.name_to_taxa[name]

    def get_taxa_str(self, idx):
        tax = self.taxa.index[idx].taxonomy
        family = tax["family"]
        genus = tax["genus"]
        species = tax["species"]

        if genus == "NA":
            return "{}**".format(family)
        elif species == "NA":
            return "{}, {}*".format(family, genus)
        else:
            return "{}, {} {}".format(family, genus, species)

    def get_taxa_str_long(self, idx):
        return "{}\n[{}]".format(self.get_taxa(idx).name, self.get_taxa_str(idx))

    def get_clustering(self):
        if self.clustering is None:
            self.clustering = self.mcmc.graph[STRNAMES.CLUSTERING_OBJ]
            for cidx, cluster in enumerate(self.clustering):
                cluster.idx = cidx
        return self.clustering

    def get_clustered_interactions(self):
        clusters = self.get_clustering()
        otu_interactions = self.get_interactions()
        cluster_interactions = np.zeros(
            shape=(
                otu_interactions.shape[0],
                len(clusters),
                len(clusters)
            ),
            dtype=float
        )
        cluster_reps = [
            next(iter(cluster.members)) for cluster in clusters
        ]
        for i in range(cluster_interactions.shape[0]):
            cluster_interactions[i] = otu_interactions[i][np.ix_(cluster_reps, cluster_reps)]
        return cluster_interactions

class KeystonenessFigure():
    def __init__(self, dataset_name: str, mcmc_pickle_path, subjset_path, fwsim_path):
        print("Loading pickle files.")
        self.md = MdsineOutput(dataset_name, mcmc_pickle_path)
        self.study = md2.Study.load(subjset_path)

        print("Loading dataframe from disk.")
        self.fwsim_df = pd.read_hdf(fwsim_path, key='df', mode='r')

        print("Compiling dataframe.")
        self.ky_df = self.generate_keystoneness_df()

        print("Compiling abundance data.")
        self.abundance_array = self.get_abundance_array()

        print("Extracting keystoneness values.")
        self.ky_array = self.get_ky_array()

        print("Extracting baseline abundances.")
        self.day20_array = self.get_day20_abundances()

    def generate_keystoneness_df(self):
        md = self.md
        fwsim_df = self.fwsim_df

        nonmembers_df = cluster_nonmembership_df(md)

        baseline = fwsim_df.loc[fwsim_df["ExcludedCluster"] == "None"]

        altered = fwsim_df.loc[fwsim_df["ExcludedCluster"] != "None"]
        altered = altered.merge(
            right=nonmembers_df,
            how="inner",
            left_on=["ExcludedCluster", "OTU"],
            right_on=["ClusterID", "OTU"]
        )

        merged = altered.merge(
            baseline[["OTU", "SampleIdx", "StableState"]],
            how="left",
            left_on=["OTU", "SampleIdx"],
            right_on=["OTU", "SampleIdx"],
            suffixes=["", "Base"]
        )

        merged["DiffStableState"] = np.log10(merged["StableStateBase"] + 1e5) - np.log10(merged["StableState"] + 1e5)

        return merged[
            ["ExcludedCluster", "SampleIdx", "DiffStableState"]
        ].groupby(
            ["ExcludedCluster", "SampleIdx"]
        ).mean().rename(columns={"DiffStableState": "Ky"})

    def get_abundance_array(self):
        md = self.md
        fwsim_df = self.fwsim_df

        clustering = md.get_clustering()
        membership_df = cluster_membership_df(md)
        merged_df = fwsim_df.merge(
            membership_df,
            how="left",
            left_on="OTU",
            right_on="OTU"
        )

        abund_array = np.zeros(shape=(len(clustering) + 1, len(clustering)))

        # Baseline abundances (no cluster removed) -- sum across OTUs (per sample), median across samples.
        subset_df = merged_df.loc[merged_df["ExcludedCluster"] == "None"]
        subset_df = subset_df[
            ["ClusterOfOTU", "SampleIdx", "StableState"]
        ].groupby(
            ["ClusterOfOTU", "SampleIdx"]
        ).sum(
            # Aggregate over OTUs  (e.g. Baseline abundance of a cluster is the sum of its constituents.)
        ).groupby(
            level=0
        ).median(
            # Aggregate over samples
        )
        for cluster in clustering:
            abund_array[0, cluster.idx] = subset_df.loc[cluster.id]

        # Altered abundances (remove 1 cluster at a time)
        for removed_cluster in tqdm(clustering, total=len(clustering), desc="Heatmap Abundances"):
            subset_df = merged_df.loc[merged_df["ExcludedCluster"] == removed_cluster.id]

            # Compute the total abundance (over OTUs) for each cluster, for each sample. Then aggregate (median) across samples.
            subset_df = subset_df[
                ["ClusterOfOTU", "SampleIdx", "StableState"]
            ].groupby(
                ["ClusterOfOTU", "SampleIdx"]
            ).sum(
                # Aggregate over OTUs
            ).groupby(
                level=0
            ).median(
                # Aggregate over samples
            )

            for cluster in clustering:
                abund_array[removed_cluster.idx + 1, cluster.idx] = subset_df.loc[cluster.id]
        return abund_array

    def get_ky_array(self):
        # Group by Cluster, aggregate (mean/median) across samples.
        agg_ky_df = self.ky_df.groupby(level=0).median()
        return np.array(
            [agg_ky_df.loc[cluster.id, "Ky"] for cluster in self.md.get_clustering()]
        )

    def get_day20_abundances(self):
        M = self.study.matrix(dtype='abs', agg='mean', times='intersection', qpcr_unnormalize=True)
        day20_state = M[:, 19]
        cluster_day20_abundances = np.zeros(len(self.md.get_clustering()))

        for cidx, cluster in enumerate(self.md.get_clustering()):
            cluster_day20_abundances[cidx] = np.sum(
                [day20_state[oidx] for oidx in cluster.members]
            )
        return cluster_day20_abundances

    def plot(self, fig):
        md = self.md
        abund_array = self.abundance_array
        ky_array = self.ky_array
        day20_array = self.day20_array

        # Main abundance grid shows the _difference_ from baseline, instead of the abundances itself.
        n_clusters = len(ky_array)

        # =========== Pre-sorting. ===========
        ky_order = np.argsort(ky_array)
        ky_order = ky_order[::-1]
        ky_array = ky_array[ky_order]

        day20_array = day20_array[ky_order].reshape(1, len(day20_array))

        baseline_array = abund_array[[0], :]
        baseline_array = baseline_array[:, ky_order]

        altered_array = abund_array[1 + ky_order, :]  # Reorder the rows first (exclude the baseline row),
        altered_array = altered_array[:, ky_order]  # Then reorder the columns.

        baseline_diff_array = np.log10(baseline_array + 1e5) - np.log10(altered_array + 1e5)
        for i in range(baseline_diff_array.shape[0]):
            baseline_diff_array[i, i] = np.nan

        # =========== Heatmap settings. ========
        gridspec_kw = {"height_ratios":[1, 1, abund_array.shape[0] - 1], "width_ratios" : [1, abund_array.shape[1]]}

        # Colors and normalization (abund)
        abund_min = np.max([
            np.min(abund_array[abund_array > 0]),
            1e5
        ])
        abund_max = np.min([
            np.max(abund_array[abund_array > 0]),
            1e13
        ])
        print("abund_min = {}, abund_max = {}".format(abund_min, abund_max))

        abund_cmap = create_cmap("Greens", nan_value="white")
        abund_norm = matplotlib.colors.LogNorm(vmin=abund_min, vmax=abund_max)

        # Colors and normalization (Ky)
        gray = np.array([0.95, 0.95, 0.95, 1.0])
        red = np.array([1.0, 0.0, 0.0, 1.0])
        blue = np.array([0.0, 0.0, 1.0, 1.0])
        n_interp=128

        top = blue
        bottom = red
        top_middle = 0.05 * top + 0.95 * gray
        bottom_middle = 0.05 * bottom + 0.95 * gray

        ky_cmap = matplotlib.colors.ListedColormap(
            np.vstack(
                [(1-t) * bottom + t * bottom_middle for t in np.linspace(0, 1, n_interp)]
                +
                [(1-t) * top_middle + t * top for t in np.linspace(0, 1, n_interp)]
            ),
            name='Keystoneness'
        )
        ky_cmap.set_bad(color="white")


        ky_min = 0.90 * np.min(ky_array) + 0.10 * np.min(baseline_diff_array[altered_array > 0])
        ky_max = 0.90 * np.max(ky_array) + 0.10 * np.max(baseline_diff_array[altered_array > 0])

        ky_norm = matplotlib.colors.TwoSlopeNorm(vmin=ky_min, vcenter=0, vmax=ky_max)
        def _forward(x):
            y = x.copy()
            positive_part = x[x > 0]
            y[x > 0] = np.sqrt(positive_part / ky_max)

            negative_part = x[x < 0]
            y[x < 0] = -np.sqrt(np.abs(negative_part / ky_min))
            return y
        def _reverse(x):
            y = x.copy()
            positive_part = x[x > 0]
            y[x > 0] = ky_max * np.power(positive_part, 2)

            negative_part = x[x < 0]
            y[x < 0] = -np.abs(ky_min) * np.power(negative_part, 2)
            return y
        ky_norm = mcolors.FuncNorm((_forward, _reverse), vmin=ky_min, vmax=ky_max)

        # Seaborn Heatmap Kwargs
        abund_heatmapkws = dict(square=False,
                                cbar=False,
                                cmap=abund_cmap,
                                linewidths=0.5,
                                norm=abund_norm)
        ky_heatmapkws = dict(square=False, cbar=False, cmap=ky_cmap, linewidths=0.5, norm=ky_norm)

        # ========== Plot layout ===========
    #     [left, bottom, width, height]
        main_x = 0.67
        main_y = 0.5
        box_unit = 0.03
        main_width = box_unit * n_clusters
        main_height = main_width
        main_left = main_x - 0.5 * main_width
        main_bottom = main_y - 0.5 * main_width
        print("Left: {}, bottom: {}, width: {}, height: {}".format(main_left, main_bottom, main_width, main_height))
        print("Right: {}, Top: {}".format(main_left + main_width, main_bottom + main_height))

        ky_ax = fig.add_axes([main_left + main_width + 0.5 * box_unit, main_bottom, box_unit, main_height])
        abundances_ax = fig.add_axes([main_left, main_bottom, main_width, main_height])
        obs_ax = fig.add_axes([main_left, main_bottom + main_height + 1.5 * box_unit, box_unit * n_clusters, box_unit])
        baseline_ax = fig.add_axes([main_left, main_bottom + main_height + 0.5 * box_unit, box_unit * n_clusters, box_unit])

        # ========= Rendering. ==========
        # ====== Bottom left: Keystoneness
        hmap_ky = sns.heatmap(
            ky_array.reshape(len(ky_array), 1),
            ax=ky_ax,
            xticklabels=False,
            yticklabels=False,
            **ky_heatmapkws
        )
        hmap_ky.xaxis.set_tick_params(width=0)
        fig.text(main_left + main_width + 2*box_unit, main_y, "Keystoneness", ha='center', va='center', rotation=-90,
                fontsize="xx-large", fontweight="bold")

        for _, spine in hmap_ky.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(1.0)

        # ====== Top right 1: Observed levels (day 20)
        hmap_day20_abund = sns.heatmap(day20_array,
                                       ax=obs_ax,
                                       xticklabels=False,
                                       yticklabels=["Observation"],
                                       **abund_heatmapkws)
        hmap_day20_abund.set_yticklabels(hmap_day20_abund.get_yticklabels(), rotation=0, fontsize="x-large")
        for _, spine in hmap_day20_abund.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(1.0)
        fig.text(main_x, main_bottom + main_height + 3 * box_unit, "Steady State Abundance", ha='center', va='center',
                fontsize="xx-large", fontweight="bold")

        # ====== Top right 2: Baseline abundances
        hmap_base_abund = sns.heatmap(baseline_array,
                                      ax=baseline_ax,
                                      xticklabels=False,
                                      yticklabels=["Simulation"],
                                      **abund_heatmapkws)
        hmap_base_abund.set_yticklabels(hmap_base_abund.get_yticklabels(), rotation=0, fontsize="x-large")
        for _, spine in hmap_base_abund.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(1.0)

        # ====== Bottom right: Abundances with clusters removed.
        ticklabels = [
            "{}{}".format(
                "H" if md.dataset_name == "Healthy" else "D",
                c_idx + 1
            )
            for c_idx in ky_order
        ]
        hmap_removed_cluster_abund = sns.heatmap(
            baseline_diff_array,
            ax=abundances_ax,
            xticklabels=ticklabels,
            yticklabels=ticklabels,
            **ky_heatmapkws
        )
        # Draw a marker ("X") on top of NaNs.
        abundances_ax.scatter(*np.argwhere(np.isnan(baseline_diff_array.T)).T + 0.5, marker="x", color="black", s=100)
        abundances_ax.set_ylabel("Module Removed", fontsize="xx-large", fontweight="bold")
        abundances_ax.set_xlabel("Per-Module Change", fontsize="xx-large", fontweight="bold")
        for _, spine in hmap_removed_cluster_abund.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(1.0)
        hmap_removed_cluster_abund.xaxis.set_ticks_position('bottom')
        hmap_removed_cluster_abund.set_xticklabels(
            hmap_removed_cluster_abund.get_xticklabels(), rotation=90, horizontalalignment='center',
            fontsize="x-large"
        )
        hmap_removed_cluster_abund.set_yticklabels(hmap_removed_cluster_abund.get_xticklabels(),
            rotation=0, fontsize="x-large")
        abundances_ax.tick_params(direction='out', length=0, width=0)

        # ======= Draw the colormaps ========
        cbar_from_main = 0.2
        cbar_width = 0.01
        cbar_height = 0.35

        # Cbar on the right (steady state diff, green)
        cax = fig.add_axes([main_left - cbar_from_main, main_y - 0.5 * cbar_height, cbar_width, cbar_height])
        sm = matplotlib.cm.ScalarMappable(cmap=abund_cmap, norm=abund_norm)
        sm.set_array(np.array([]))
        cbar = fig.colorbar(sm, cax=cax)

        yticks = cbar.get_ticks()
        yticklabels = [str(np.log10(y)) for y in yticks]
        yticklabels[0] = "<{}".format(yticklabels[0])
        cax.set_yticklabels(yticklabels, fontsize="x-large")
        cax.set_ylabel("Log-Abundance", fontsize="xx-large", fontweight="bold")

        # Cbar on the left (Keyst., RdBu)
        cax = fig.add_axes([main_left - cbar_from_main - 2*cbar_width, main_y - 0.5 * cbar_height, cbar_width, cbar_height])
        sm = matplotlib.cm.ScalarMappable(cmap=ky_cmap, norm=ky_norm)
        sm.set_array(np.array([]))
        cbar = fig.colorbar(sm, cax=cax)
        cax.yaxis.set_ticks_position('left')
        cax.set_ylabel("Log-Difference from Base", fontsize="xx-large", fontweight="bold")
        cax.yaxis.set_label_position("left")

        yticks = cbar.get_ticks()
        yticklabels = ["{:.1f}".format(y) for y in yticks]
        cax.set_yticklabels(yticklabels, fontsize="x-large")

        # Legend label text
        fig.text(
            main_left - cbar_from_main - cbar_width,
            main_y + 0.5 * cbar_height + 0.05,
            "Legend",
            ha='center', va='center',
            fontweight='bold', fontsize="xx-large"
        )
