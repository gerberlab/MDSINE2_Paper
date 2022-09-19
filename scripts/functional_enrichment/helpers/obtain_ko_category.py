from Bio.KEGG import REST
from maps import *
from pathlib import Path

import argparse


def extract_all_list(database_name):
    """
    extract all entries in a given database
    :param (str) database_name: name of a database
    return : (dict) (str) entry code -> (str) entry name
    """

    all_entries = REST.kegg_list(database_name).read()

    return all_entries.strip()


def hash_entries(entries):
    """
    has the contents of the given string
    :param (str) entries: the entire entries of a database
    """
    entries_li = entries.split("\n")
    hash_table = {}
    for entry in entries_li:
        entry_split = entry.split("\t")
        key = entry_split[0].split(":")[-1]
        value = entry_split[1].split(",")[0]
        hash_table[key] = value
    return hash_table


def get_gene_organism_info():
    """
    return : (dict) (str) gene name -> (str) organism type
    """

    all_genes = REST.kegg_list("organism").read().strip().split("\n")[0:-1]

    gene_organism_dict = {}
    for gene in all_genes:
        gene_info_split = gene.split("\t")
        gene_name = gene_info_split[1]
        organism = gene_info_split[-1].split(";")[0]
        gene_organism_dict[gene_name] = organism.lower()

    return gene_organism_dict


def classify_ko_organism(ko_list, gene_organism_dict):
    """
    classify whether a KO is associated with a prokaryote or eukaryote or both organisms
    """
    pass


def get_ko_gene_information(ko_name):
    """
    :param (str) ko_name: Kegg Orthology
    return (List[str]) genes associated with a KO
    """

    ko_info = REST.kegg_get(ko_name).read().split("\n")
    i = 0
    found_gene = False
    genes = []
    for w in ko_info:
        if "GENES" in w:
            found_gene = True
            gene_idx = i
            genes.append(w.split()[1].lower().strip(":"))
        else:
            if found_gene:
                if len(w[0:3].strip()) == 0:
                    genes.append(w.split()[0].lower().strip(":"))
                else:
                    found_gene = False
    return genes


def get_ko_category(ko, gene_category_dict):
    """
    classify KO as prokaryotic or eukaryotic or both
    :param (List[str]) ko_list: names of KO
    :param (dict) gene_category_dict: (str) gene name ->(str) prokaryotic / eukaryotic

    return (dict): (str) ko name -> (str)  prokaryotic / eukaryotic / both
    """

    is_prokar = False
    is_eukar = False
    ko_genes = get_ko_gene_information(ko)
    try:
        for gene in ko_genes:
            gene_cat = gene_category_dict[gene]
            if gene_cat == "prokaryotes":
                is_prokar = True
            elif gene_cat == "eukaryotes":
                is_eukar = True
            else:
                print("different type, gene: {}, ko: {}".format(
                    gene, ko))
            if is_prokar and is_eukar:
                return "both"
    except Exception:
        print("issue with gene: {}, KO:{}".format(gene, ko))

    if is_prokar:
        return "prokaryote"
    else:
        return "eukaryote"


def obtain_functional_unit_category(associated_kos, ko_info_dict):
    """
    :param (List[str]) : names of ko's associated with a functional unit
    :param (dict) ko_info_dict: (str) ko -> (str) prokaryote or eukaryote or both

    return: (dict) (str) category type -> (int) count for each category
    """

    dict_ = {"prokaryote": 0, "eukaryote": 0, "both": 0}
    for ko in associated_kos:
        category = ko_info_dict[ko]
        if category not in dict_:
            raise Exception("Wrong cateogry")
        else:
            dict_[category] += 1

    return dict_


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--funit_ko_file", help="the text file containing information"
                                                "about functional units and their associated KOs")
    parser.add_argument("--save_loc", help="location where the dictionary is saved "
                                           "as a pickle file")
    parser.add_argument("--save_name", help="name of the pickle file")

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_arguments()
    all_KOs = extract_all_list("orthology")
    ko_details_dict = hash_entries(all_KOs)
    gene_organisms_dict = get_gene_organism_info()
    print(gene_organisms_dict)

    #manual fix
    gene_organisms_dict["ag"] = "prokaryotes"

    ko_category_dict = {}
    i = 0
    for ko in ko_details_dict:
        ko_category = get_ko_category(ko, gene_organisms_dict)
        ko_category_dict[ko] = ko_category
        i += 1
        if i % 200 == 0:
            print(i)

    df_series = pd.Series(ko_category_dict)
    save_loc = Path(args.save_loc)
    save_loc.mkdir(parents=True, exist_ok=True)

    df_series.to_csv(save_loc / "ko_category.csv")
    """
    funit_trait_path = Path(args.funit_trait_path) if ".txt" in args.funit_trait_path else \
        Path(args.funit_trait_path) + ".txt"
    funit_trait_path.mkdir(parents=True, exist_ok=True)
    dict_funit_ko = map_functional_units_to_traits(funit_trait_path)

    funit_category_count_dict = {}
    for funit in dict_funit_ko:
        funit_category_count_dict[funit] = obtain_functional_unit_category(
            dict_funit_ko[funit], ko_category_dict)
        print(funit, funit_category_count_dict[funit] )
    """



