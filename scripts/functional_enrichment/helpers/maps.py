import pandas as pd
import numpy as np


def map_traits_to_sequence_units(table_loc, bugs_of_interest, sep=","):
    """
    establish a link between traits and associated sequence units
    :param (Path) table_loc: location of the table illustrating the predicted abudance for
                             each sequence variant
    :param (dict) bugs_of_interest : (str) bug_name -> (int) index of the bug
    :param (str) sep:
    returns: (dict) (str) trait_name -> (List[str]) names of sequence units (OTU or ASV)
    """

    table = pd.read_csv(table_loc, sep=sep, index_col=0)
    table_numpy = table.to_numpy()
    bug_names = table.index
    trait_names = table.columns
    trait_sequence_map = {}
    for i in range(len(trait_names)):
        trait = trait_names[i]
        for j in range(len(bug_names)):
            if bug_names[j] in bugs_of_interest and table_numpy[j][i] != 0:
                if trait not in trait_sequence_map:
                    trait_sequence_map[trait] = []
                trait_sequence_map[trait].append(bug_names[j])

    return trait_sequence_map


def update_map(dict_, key, value):
    """
    updates the dictionary by adding the value to the list associated with key
    :param (dict) dict_: (str/int) -> List[str/int]
    :param (str/int) key: key in the dictionary
    :param (str/int) value: value to tbe added in the dictionary
    """
    if key not in dict_:
        dict_[key] = [value]
    else:
        dict_[key].append(value)


def run_filter0(funit_object):
    """
    filter to select modules that are dominantly prokaryotic
    :param (PathwayModule) funit_object: instantiation of the PathwayModule object
    returns: (bool) whether or not the module/pathway passes filter 0
    """
    if funit_object.get_eukaryote_count() > funit_object.get_prokaryote_count() :
        return False
    else:
        return True


def map_functional_units_to_traits(txt_file_loc, valid_trait_dict, funit_info_dict=None):
    """
    establish a link between functional units and the traits associated with the units
    :param (Path) txt_file_loc: path to the txt file
    :param (dict) funit_info_dict: (str) pathway/module name -> (Module) Module object containing
                                              relevant information
    :param (dict) valid_trait_dict: (str) name of trait -> (object)
    returns (dict): (str) functional unit name -> (List[str]) associated traits
    """

    enable_filter0 = False # need to add filters
    observed = {}
    repeated = {}
    funit_trait_map = {}
    file = open(txt_file_loc, "r")
    for line in file:
        line_split = line.strip().split()
        if len(line_split) > 1:
            group = line_split[0]
            if group not in observed:
                observed[group] = 1
            else:
                if group not in repeated:
                    repeated[group] = 0
                repeated[group] += 1
            if group in valid_trait_dict:
                for i in range(1, len(line_split)):
                    f_unit = line_split[i]
                    if funit_info_dict is not None:
                        if run_filter0(funit_info_dict[f_unit]):
                            update_map(funit_trait_map, f_unit, group)
                    else:
                        update_map(funit_trait_map, f_unit, group)
    return funit_trait_map


def map_two_dicts(dict1, dict2):
    """
    link the key of one dict with value of another; the values in dict2 must
    be the key of dict2
    :param (dict) dict1: (str/int) -> List[str/int]
    :param (dict) dict2: (str/int) -> List[str/int]

    return: (dict) (str/int) -> List[str/int]
    """

    dict_ = {}
    for key1 in dict1:
        dict_[key1] = set()
        for value1 in dict1[key1]:
            if value1 in dict2:
                dict_[key1] = dict_[key1].union(set(dict2[value1]))

    final_dict = {}
    for key in dict_:
        final_dict[key] = list(dict_[key])

    return final_dict


def export_results(enrich_results, non_corrected_p_values, funits, funit_names, loc,
                   f_name):
    """
    save the results in  a txt file

    :param (tuple([bool], [float])) enrich_results: enrichment results and corrected p_vals
    :param (List[float]) non_corrected_p_values:  non-corrected p_values
    :param (List[str]) funits: names functional units
    :param (dict) funit_names: (str) functional unit short name -> (str) full name of functional units
    :param (Path) loc: location where the file is saved
    :param (str) f_name: name of the file that is exported
    """

    enriched_funits = []
    enriched_p_values = []
    enriched_adjusted_p = []

    for j in range(len(enrich_results[0])):
        if enrich_results[0][j]:
            enriched_funits.append(funits[j])
            enriched_p_values.append(non_corrected_p_values[j])
            enriched_adjusted_p.append(enrich_results[1][j])

    enriched_cluster_module = {}
    enriched_p = {}
    enriched_adj_p = {}

    for i in range(len(enriched_funits)):
        info = enriched_funits[i].split()
        cluster_id = info[1]
        if cluster_id not in enriched_cluster_module:
            enriched_cluster_module[cluster_id] = []
        if cluster_id not in enriched_p:
            enriched_p[cluster_id] = []
        if cluster_id not in enriched_adj_p:
            enriched_adj_p[cluster_id] = []
        name = info[2]
        enriched_cluster_module[cluster_id].append(funit_names[name])
        enriched_p[cluster_id].append(enriched_p_values[i])
        enriched_adj_p[cluster_id].append(enriched_adjusted_p[i])

    all_lines = "Total Enrichment: " + str(len(enriched_funits)) + "\n"
    for keys in enriched_cluster_module:
        all_lines = all_lines + "Cluster " + keys + "\n"
        all_lines = all_lines + "module name" + "\t" + "p value" + "\t" + "adjusted_p" + "\n"
        n = len(enriched_cluster_module[keys])
        for i in range(n):
            name = enriched_cluster_module[keys][i]
            p1 = enriched_p[keys][i]
            p2 = enriched_adj_p[keys][i]
            all_lines += name + "\t" + str(p1) + "\t" + str(p2) + "\n"
        all_lines = all_lines + "\n"
    results = open(loc/"{}_results.txt".format(f_name), "w")
    results.write(all_lines)