from Bio.KEGG import REST
from KO_scraper import KO, extract_all_list, hash_entries
from pathlib import Path

import argparse
import pickle
import pandas as pd


class PathwayModule:
    """
    base class for Kegg Pathway / Module
    Attributes:
        (str) name: name of the pathway/module
        (str) detailed name: detailed name of the module / pathway
        (List[str]) associated_kos: KOs related to the pathway/module
    """

    def __init__(self, name):
        self.name = name
        self.detailed_name = ""
        self.associated_kos = []
        self.prokaryotes_count = 0
        self.eukaryotes_count = 0
        self.both_count = 0

    def add_ko(self, ko):
        self.associated_kos.append(ko)

    def get_associated_kos(self):
        return self.associated_kos

    def increase_count(self, org_type):
        if org_type == "prokaryote":
            self.prokaryotes_count += 1
        elif org_type == "eukaryote":
            self.eukaryotes_count += 1
        elif org_type == "both":
            self.both_count += 1
        else:
            print("{} not a valid organism type".format(org_type))

    def get_prokaryote_count(self):
        return self.prokaryotes_count

    def get_eukaryote_count(self):
        return self.eukaryotes_count

    def get_both_count(self):
        return self.both_count


def export_csv(hash_table, loc, name):
    """
    export the contents of hash tables to a csv file
    :param (dict) hash_table: (str) symbol -> (str) name
    :param (Path) loc: location where the csv file is exported
    :param (str) name: name of the csv file
    """

    df = pd.Series(hash_table)
    df.to_csv(loc/"{}.csv".format(name), header=False)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ko_details_dir", help="directory of a  dictionary whose key "
                        "is a KO name and value is the corresponding KO object. The dictionary"
                                                 "is store as a pkl file")
    parser.add_argument("--save_loc", help='location where the output files are saved')
    parser.add_argument("--kegg_type", help="pathway or module, the KEGG entry")

    return parser.parse_args()


def export_count_info(module_info, loc, name):
    """
    export information about the count of organism associated with the module as a txt file;
    useful for debugging purposes

    :param (dict) module_info: (str) module_name -> (PathwayModule) PathwayModule object
    :param (Path) loc: location where the output is saved
    :param (str) name: name of the output file
    """
    all_str = ""
    for module in module_info:
        info = module_info[module]
        all_str = all_str + "{}, eu:{}, pro:{}, both:{}\n".format(module, info.get_eukaryote_count(),
                            info.get_prokaryote_count(), info.get_both_count())

    file = open(loc / "{}.txt".format(name), "w")
    file.write(all_str)
    file.close()


if __name__ == "__main__":
    args = parse_arguments()
    ko_details_dict = pickle.load(open(args.ko_details_dir, "rb"))
    all_entries = extract_all_list(args.kegg_type)
    hash_table_entries = hash_entries(all_entries.strip())

    save_loc = Path(args.save_loc)
    save_loc.mkdir(parents=True, exist_ok=True)
    export_csv(hash_table_entries, save_loc, args.kegg_type + "s_names")

    all_pathway_modules = {}
    for ko in ko_details_dict:
        ko_obj = ko_details_dict[ko]
        funit_ko = ko_obj.get_associated_modules() if args.kegg_type == "module" else \
            ko_obj.get_associated_pathways()
        for funit in funit_ko:
            if funit not in all_pathway_modules:
                all_pathway_modules[funit] = PathwayModule(funit)
            all_pathway_modules[funit].add_ko(ko)
            all_pathway_modules[funit].increase_count(ko_obj.get_ko_category())
    export_count_info(all_pathway_modules, save_loc, "{}_category_counts".format(args.kegg_type))