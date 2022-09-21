from maps import *
import argparse
import mdsine2 as md2
import copy
import pickle

from pathlib import Path
from mdsine2.names import STRNAMES
from scipy.special import comb
from scipy.stats import hypergeom
from statsmodels.stats import multitest as mtest
from pathway_module_scraper import PathwayModule


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait_abundance_table_path", help="table showing the abundance of KOs for each "
                        "sequence units", dest="trait_table_path")
    parser.add_argument("--funit_trait_txt_file_path", help="text file that includes information about "
                        "the relationship between trait groups and functional units",
                        dest="funit_trait_path")
    parser.add_argument("--mcmc_pkl_loc", dest="mcmc_loc", help="location of the mcmc simulations")
    parser.add_argument("--filter1", type=float, help="percentage threshold")
    parser.add_argument("--filter2", type=int, help = "min number of bugs in a cluster that must be annotated\
             to a functional unit of interested ")
    parser.add_argument("-f3", "--filter3", type = int, required=True,
                        help="minimum number of KOs/ECs that must be associated with a functional\
                        unit to be considered valid for enrichment")
    parser.add_argument("--names_file_loc", help="name of the folder where the files "
                        "containing the names of functional units are located")
    parser.add_argument("--pathway_module_pkl_file", default=None, help="location of the pickle file"
                        "containing details about pathway or module")
    parser.add_argument("--funit_type", help="pathway/module/Cazyme, the functional unit being"
                                             "tested for enrichment")
    parser.add_argument("--save_loc", help="location where the data frames are saved")
    return parser.parse_args()


def extract_cluster_information(mcmc):
    """
    obtain information abouts clusters from mcmc pkl file

    :param (BaseMCMC) mcmc: results of the mcmc inference stored as a pkl file
    return
    """
    def clusterize(labels, taxa_list):
        """returns the consenus cluster as a dictionary, (int) cluster_id ->
        (List[str]) ids of OTUs belonging to the cluster

        :param labels: (List[int]) the consensus cluster id of bugs in taxa_list
        :param taxa_list: (List[str]) ids of OTU
        """
        cluster = {}
        for i in range(len(labels)):
            if labels[i] + 1 not in cluster:
                cluster[labels[i] + 1] = []
            cluster[labels[i] + 1].append(taxa_list[i].replace("OTU", "ASV"))
        return cluster

    clustering = mcmc.graph[STRNAMES.CLUSTERING].clustering
    consensus_cluster_labels = clustering.toarray()
    taxa_list = []
    taxas = mcmc.graph.data.taxa
    for taxa in taxas:
        taxa_list.append(taxa.name)
    cluster_dict = clusterize(consensus_cluster_labels, taxa_list)

    return cluster_dict


def compute_p_value(N, M, n, k):
    """
    computes the hypergeometric p-values
    :params (int) N, M, n, K: parameters for hyper-geometric distribution
    return: (float) p-value
    """
    lim = min(n, M)
    p = 0
    for i in range(k, lim + 1):
        p += comb(M, i) * comb(N - M, n - i) / comb(N, n)

    #p1 = hypergeom.sf(k-1, N, M, n) #to check the accuracy of p

    return p


def perform_enrichment_analysis(cluster_dict, funit_bug_dict, funit_trait_dict,
                                percent_threshold, min_bugs, min_trait_size):
    """
    performs enrichment analysis for each functional unit for each cluster
    :param (dict) cluster_dict: (str) cluster id -> List[str] list of OTU names
    :param funit_bug_dict: (str) functional unit name -> List[str] list of OTU names
    :param funit_trait_dict: (str) functional unit names -> List[str] list of traits
    :param (float) percent_threshold: upper limit on the percentage of bugs associated with a module
    :param (int) min_bugs: min number of annotated bugs in a cluster
    :param (int) min_trait_size: min number of traits (EC/KO) related to a functional unit
    """

    cluster_enriched_p = {}  # p values for functional units tested for enrichment
    cluster_enriched_funit = {}  # functional units tested for enrichment in a cluster
    cluster_all_p = {}  # cluster wise p values for all functional units
    cluster_all_funit = {}  # all functional units associated with a cluster
    total_bugs = sum([len(cluster_dict[id_]) for id_ in cluster_dict])
    max_n_funit = percent_threshold * total_bugs

    for id in cluster_dict:
        n_bugs_cluster = len(cluster_dict[id])
        cluster_enriched_p[id] = []
        cluster_enriched_funit[id] = []
        cluster_all_p[id] = []
        cluster_all_funit[id] = []
        for funit in funit_bug_dict:
            if len(funit_bug_dict[funit]) != 0:
                if len(funit_bug_dict[funit]) < max_n_funit and len(funit_trait_dict[funit]) > min_trait_size:
                    n_bugs_annotated = len(funit_bug_dict[funit])
                    n_bugs_cluster_annotated = 0
                    for bug in cluster_dict[id]:
                        if bug in funit_bug_dict[funit]:
                            n_bugs_cluster_annotated += 1
                    p_val = compute_p_value(total_bugs, n_bugs_annotated,
                                           n_bugs_cluster, n_bugs_cluster_annotated)
                    cluster_all_p[id].append(p_val)
                    cluster_all_funit[id].append(funit)
                    if n_bugs_cluster_annotated > min_bugs:
                        cluster_enriched_p[id].append(p_val)
                        cluster_enriched_funit[id].append(funit)
        return cluster_enriched_p, cluster_enriched_funit, cluster_all_p, cluster_all_funit


def get_names(csv_file):

    """
    get the names of the pathway or module or cazyme
    :param (Path) csv_file: location + name pf the csv file

    return: (dict) (str) KM/KP/Cazy code -> (str) name
    """

    name_file = pd.read_csv(csv_file, sep=",", header=None).values
    name_dict = {}
    for row in name_file:
        name_dict[row[0]] = row[0] + ", " + row[1].split(",")[0]

    return name_dict


def pivot_df(values_dict, funit_dict, enriched_res, names, loc, f_name):
    """
    arrange the data in format that's best for plotting and save it
    :param (dict) values_dict:  (str) cluster_id -> (List[float]) p-values
    :param  (dict) funit_dict: (str) cluster_id to (List[str]) KM/KP
    :param  (List[List[bool], List[float]) enriched_res: results of significance test and p-values
    :param  (dict) names: (str) KM / KP -> (str) corresponding name
    :param  (str) f_name: name of the pickle file

    returns: (pd.DataFrame)
    """

    cluster_row = []
    module_names = []
    values = []
    count = 0

    for keys in funit_dict:
        for i in range(len(funit_dict[keys])):
            if enriched_res[1][count] < 0.05:
                module_names.append(names[funit_dict[keys][i]])
                values.append(enriched_res[1][count])
                cluster_row.append("Cluster " + str(keys))
            count += 1

    data_frame = pd.DataFrame({"cluster_id": cluster_row, "module_name": module_names,
                               "p_value": values})
    df_pivot = data_frame.pivot(values="p_value", index="module_name",
                                columns="cluster_id")
    df_pivot = df_pivot.fillna(1)
    df_pivot.to_pickle(loc / "{}.pkl".format(f_name))

    return df_pivot


if __name__ == "__main__":

    args = parse_arguments()
    mcmc = md2.BaseMCMC.load(args.mcmc_loc)
    taxas = mcmc.graph.data.taxa
    bugs_mcmc_idx = {taxas[i].name.replace("OTU", "ASV") : i  for i in range(len(taxas))}
    dict_trait_asv = map_traits_to_sequence_units(Path(args.trait_table_path),
                                                  bugs_mcmc_idx, "\t")
    pathway_module_info=None
    if args.pathway_module_pkl_file is not None:
        pathway_module_info = pickle.load(open(args.pathway_module_pkl_file, "rb"))

    funit_trait_path = Path(args.funit_trait_path) if ".txt" in args.funit_trait_path else \
        Path(args.funit_trait_path) + ".txt"
    dict_funit_trait = map_functional_units_to_traits(funit_trait_path, dict_trait_asv,
                                                      funit_info_dict=pathway_module_info)
    dict_funit_asv = map_two_dicts(dict_funit_trait, dict_trait_asv)
    cluster = extract_cluster_information(mcmc)

    valid_p, valid_funits, all_p, all_funits = perform_enrichment_analysis(cluster, dict_funit_asv,
                                                                           dict_funit_trait, args.filter1,
                                                                           args.filter2, args.filter3)
    all_valid_p = [] #gather all p values in a single list
    all_valid_funits = [] #corresponding functional units
    for keys in valid_p:
        all_valid_p = all_valid_p + valid_p[keys]
        modules_cluster = []
        for k in valid_funits[keys]:
            modules_cluster.append("Cluster " + str(keys) + " " + k)
        all_valid_funits = all_valid_funits + modules_cluster

    adjusted_p = []

    if len(all_valid_p) != 0:
        adjusted_p = mtest.multipletests(copy.deepcopy(all_valid_p), alpha=0.05,
                                         method="fdr_bh", is_sorted=False)
    else:
        print("corrected p is of size 0")

    funit_names = get_names(args.names_file_loc)
    path = Path(args.save_loc)
    path.mkdir(exist_ok=True, parents=True)
    pivoted_df = pivot_df(valid_p, valid_funits, adjusted_p, funit_names, path,
                          "df_{}".format(args.funit_type))
    print(pivoted_df)
    export_results(adjusted_p, all_valid_p, all_valid_funits, funit_names, path,
                   "kegg_{}".format(args.funit_type))
    with open(path/"enrichment_results_{}".format(args.funit_type), "wb") as f:
        pickle.dump(pivoted_df, f)
