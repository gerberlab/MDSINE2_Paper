#creates a table with OTU names in rows and Sample ID as columns.
#The entries in the table correspond to the biomasses of OTUs at the time
#of sample collection

import numpy as np
import pickle as pkl
import pandas as pd
import argparse
from pathlib import Path

def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--study_file", required=True)
    parser.add_argument("-o", "--output_loc", required=True)
    parser.add_argument("-on", "--output_name", required=True)

    return parser.parse_args()

def create_abundance_table(subjset, save_loc, save_name):
    """create a table illutrating the abundances of OTUs in CFUs/g and saves it
       as a text file"""

    all_ = []
    subj_names = [subj.name for subj in subjset]
    taxas = subjset.taxa
    otu_names = [taxa.name for taxa in taxas]
    for name in subj_names:
        all_.append(subjset[name].matrix()["abs"])
    stacked = np.hstack(all_)
    col_names = [i+1 for i in range(stacked.shape[1])]
    index = otu_names
    df = pd.DataFrame(stacked, columns=col_names, index=index)
    df.columns.name = "#OTU ID"

    save_loc.mkdir(parents=True, exist_ok=True)

    df.to_csv(save_loc / "{}.txt".format(save_name), sep="\t",
        index=True, header=True, index_label="#OTU ID")

if __name__ == "__main__":

    args = parse_arguments()
    filename = ""
    if ".pkl" in args.study_file:
        filename = args.study_file
    else:
        filename = args.study_file + ".pkl"
    with open(filename, "rb") as f:
        subjset = pkl.load(f)
    create_abundance_table(subjset, Path(args.output_loc), args.output_name)
    print("Created abundance table")
