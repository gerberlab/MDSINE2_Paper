import argparse
import pickle
import urllib

from Bio.KEGG import REST
from pathlib import Path


class KO:
    """
    base class for Kegg Orthologs
    Attributes:
        (str) name: name of KO used in the database
        (List[str]) genes: the genes associated with the KO
        (List[str]) modules: the modules associated with the KO
        List[str]) pathways: the pathways associated with the KO
        (str) organism category: whether the KO is related to prokaryotes or eukaryotes or both
    """
    def __init__(self, name):
        self.name = name
        self.genes = []
        self.modules = []
        self.pathways = []
        self.organism_category = None

    def update_attributes_info(self, attributes_li):
        """
        :param (List[str]) attributes_li: the attributes we are interested in updating
        """

        attributes_dict = ko_parser(self.name, pathway=" ", module=" ", genes=":", disease=" ",
                                    name=None)
        for att in attributes_li:
            try:
                if att == "genes":
                    self.genes = attributes_dict["genes"]
                elif att == "module":
                    self.modules = attributes_dict["module"]
                elif att == "pathway":
                    self.pathways = attributes_dict["pathway"]
                else:
                    print("Attribute: {} not found".format(att))
            except KeyError:
                print("{} not present for {}".format(att, self.name))

    def update_ko_category(self, gene_category_dict):
        self.organism_category = find_ko_category(self.genes, self.name, gene_category_dict)

    def get_associated_genes(self):
        return self.genes

    def get_associated_pathways(self):
        return self.pathways

    def get_associated_modules(self):
        return self.modules

    def get_ko_category(self):
        return self.organism_category


def get_gene_organism_info():
    """
    return : (dict) (str) gene name -> (str) organism type
    """
    all_genes = REST.kegg_list("organism").read().strip().split("\n")[0:-1]
    gene_organism_dict = {}
    for gene in all_genes:
        gene_info_split = gene.split("\t")
        gene_name = gene_info_split[1].upper()
        organism = gene_info_split[-1].split(";")[0]
        gene_organism_dict[gene_name] = organism.lower()
    # manual fix
    gene_organism_dict["ag"] = "prokaryotes"
    return gene_organism_dict


def find_ko_category(ko_gene_list, ko, gene_info=None):
    """
    :param (List[str]) ko_gene_list: list of genes associated with the KO
    :param (str) ko: Kegg Ortholog
    :param (dict) gene_info: (str) -> (str) organism type
    """

    is_prokar = False
    is_eukar = False
    if gene_info is None:
        gene_info = get_gene_organism_info()
    try:
        for gene in ko_gene_list:
            gene_cat = gene_info[gene]
            if gene_cat == "prokaryotes":
                is_prokar = True
            elif gene_cat == "eukaryotes":
                is_eukar = True
            else:
                print("different type, gene: {}, ko: {}".format(
                    gene, ko))
            if is_prokar and is_eukar:
                return "both"
    except KeyError:
        print("issue with gene: {}, KO:{}".format(gene, ko))

    if is_prokar:
        return "prokaryote"
    else:
        return "eukaryote"


def extract_all_list(database_name):
    """
    extract all entries in a given database
    :param (str) database_name: name of a database
    return : (dict) (str) entry code -> (str) entry name
    """

    all_entries = REST.kegg_list(database_name).read()
    return all_entries.strip()


def ko_parser(ko, **separator):
    """
    parse information about a given KO
    :param (str) ko: name of Kegg Ortholog
    :param (List[str]) attributes: the KO entry attributes we want
    :param separator: (key) attribute -> (str) symbol used to separate
                            attribute name and identifier
    """
    print("parsing information about:", ko)
    try:
        ko_info = REST.kegg_get(ko).read().rstrip().split("\n")
        ko_contents = {}
        current_attribute = None
        for line in ko_info:
            line_split = line
            attribute = line_split[:12].strip()
            if attribute != "":
                current_attribute = attribute.lower()
            if current_attribute not in ko_contents:
                ko_contents[current_attribute] = []
            if current_attribute in separator:
                if separator[current_attribute] is not None:
                    attribute_identifier = line[12:].split(
                             separator[current_attribute])[0]
                    if current_attribute == "pathway":
                        attribute_identifier = attribute_identifier.replace("map", "ko")
                    ko_contents[current_attribute].append(attribute_identifier)
            else:
                ko_contents[current_attribute].append(line[12:].strip())
        return ko_contents
    except Exception:
        print("Error with {}", ko)


def hash_entries(entries):
    """
     has the contents of the KEGG entries
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


def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument("--save_loc", help="location where the pkl files are saved")
    parser.add_argument("--output_name", help="name of the pkl files used in saving KO objects")
    parser.add_argument("--start_idx", required=False, default=0)

    return parser.parse_args()


def create_ko_pathway_module_file(ko_info, save_loc, name, type_):
    """
     create a text file that lists KOs and modules/pathways associated with the KO
    :param (dict) ko_info: (str) ko name -> (KO) KO object
    :param (Path) save_loc: location where the output is saved
    :param (str) name: name of the output file
    :param (str) type_: module or pathway
    """
    details_str = ""

    for ko in ko_info:
        ko_obj = ko_info[ko]
        funit_ko = ko_obj.get_associated_modules() if type_ == "modules" else ko_obj.get_associated_pathways()
        details_str = details_str + ko + " " + " ".join(funit_ko) + "\n"

    dest_file = open(save_loc/"{}.txt".format(name), "w")
    dest_file.write(details_str)
    dest_file.close()


if __name__ == "__main__":
    all_KOs = extract_all_list("orthology")
    ko_hashed = hash_entries(all_KOs)
    ko_details_dict = {}
    gene_type_info = get_gene_organism_info()
    i = 0
    for ko in ko_hashed:
        ko_obj = KO(ko)
        ko_obj.update_attributes_info(["pathway", "module", "genes"])
        ko_obj.update_ko_category(gene_type_info)
        i += 1
        ko_details_dict[ko] = ko_obj
        #if i == 20:
         #   break

    args = parse_arguments()
    save_path = Path(args.save_loc)
    save_path.mkdir(parents=True, exist_ok=True)
    with open(save_path/"{}.pkl".format(args.output_name), "wb") as handle:
        pickle.dump(ko_details_dict, handle)

    create_ko_pathway_module_file(ko_details_dict, save_path, "kegg_modules_ko", "modules")
    create_ko_pathway_module_file(ko_details_dict, save_path, "kegg_pathways_ko", "pathways")
