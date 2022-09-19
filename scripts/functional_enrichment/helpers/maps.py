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


def map_functional_units_to_traits(txt_file_loc):
    """
    establish a link between functional units and the traits associated with the units
    :param (Path) txt_file_loc: path to the txt file
    returns (dict): (str) functional unit name -> (List[str]) associated traits
    """

    enable_filter0 = False # need to add filters
    funit_trait_map = {}
    file = open(txt_file_loc, "r")
    for line in file:
        line_split = line.strip().split()
        if len(line_split) > 1:
            group = line_split[0]
            for i in range(1, len(line_split)):
                f_unit = line_split[i]
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
