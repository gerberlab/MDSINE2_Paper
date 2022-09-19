import pandas as pd
from Bio.KEGG import REST
from pathlib import Path

def extract_all_list(database_name):
    """
    extract all entries in a given database
    :param (str) database_name: name of a database
    return : (dict) (str) entry code -> (str) entry name
    """

    all_entries = REST.kegg_list(database_name).read()

    return all_entries


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


def export_csv(hash_table, loc, name):
    """
    export the contents of hash tables to a csv file
    :param (dict) hash_table: (str) symbol -> (str) name
    :param (Path) loc: location where the csv file is exported
    :param (str) name: name of the csv file
    """

    df = pd.Series(hash_table)
    df.to_csv(loc / "{}.csv".format(name), header=False)


if __name__ == "__main__":
    '''
    path = Path("datasets/gibson/functional_enrichment/ref_files/common")
    path.mkdir(exist_ok=True, parents=True)

    all_entries_modules = extract_all_list("module")
    hash_table_modules = hash_entries(all_entries_modules.strip())
    export_csv(hash_table_modules, path, "module_names")

    all_entries_pathway = extract_all_list("pathway")
    hash_table_pathway = hash_entries(all_entries_pathway.strip())
    export_csv(hash_table_pathway, path, "pathway_names_")
    '''

    all_entries = REST.kegg_list("organism").read()
    info = REST.kegg_get("K16879").read()
    print(info)
    print(all_entries)