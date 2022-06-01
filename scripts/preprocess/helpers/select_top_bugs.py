#obtains the N most abundant bugs from a study object and saves it as a study
#object
'''
    python retain_top_bugs.py \
        -f "/Users/microbiome/Desktop/MDSINE2_Paper/analysis/output/gibson/preprocessed/gibson_healthy_agg_taxa_filtered.pkl"\
        -o "/Users/microbiome/Desktop/top_OTUs" \
        -n 25
'''
import mdsine2 as md2
import numpy as np
import argparse
from pathlib import Path
import pandas as pd

def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--study_file", required=True,
        help="The location of the study object")
    parser.add_argument("-o", "--output_loc", required=True,
       help="location where the output is saved")
    parser.add_argument("-n", "--n_bugs", required=True, type=int,
       help="the top n most abundant bugs that are retained")

    return parser.parse_args()

def count_times(df, times, times_cnts, t2idx):
    """
    count the number of times data is collected at a given sample time point
    """

    for col in df.columns :
        if col in times:
            times_cnts[t2idx[col]] += 1

    return times_cnts

def add_unequal_col_dataframes(df, dfother, times, times_cnts, t2idx):
    """combine 2 dataframe that may have different columns headings"""

    times_cnts = count_times(dfother, times, times_cnts, t2idx)
    if df is None:
        return dfother, times_cnts

    cols_toadd_df = []
    cols_toadd_dfother = []
    for col in dfother.columns:
        if col not in df.columns:
            cols_toadd_df.append(col)
    for col in df.columns:
        if col not in dfother.columns:
            cols_toadd_dfother.append(col)

    df = pd.concat([df,
        pd.DataFrame(np.zeros(shape=(len(df.index), len(cols_toadd_df))),
            index=df.index, columns=cols_toadd_df)], axis=1)
    dfother = pd.concat([dfother,
        pd.DataFrame(np.zeros(shape=(len(dfother.index), len(cols_toadd_dfother))),
            index=dfother.index, columns=cols_toadd_dfother)], axis=1)


    return dfother.reindex(df.index) + df, times_cnts

def combine_dfs_subject(subjset, dtype):
    """
       combines the abundance of OTUs in different subjects present in subjset
       together in a DataFrame
    """

    times = []
    for subj in subjset:
        times = np.append(times, subj.times)
    times = np.sort(np.unique(times))
    t2idx = {}
    for i, t in enumerate(times):
        t2idx[t] = i
    times_cnts = np.zeros(len(times))

    df = None
    for subj in subjset:
        dfnew = subj.df()[dtype]
        df, times_cnts = add_unequal_col_dataframes(df = df, dfother = dfnew,
        times = times, times_cnts = times_cnts, t2idx = t2idx)
    df = df / times_cnts
    return df

def get_top_abundance_data(study, N):
    """
       makes a study object consisting of the N most abundant OTUs
       @parameters
       -----------------------
       subjset : mdsine2.base.Study

       @returns
       -----------------------
       mdsine2.base.study
    """

    combined_df = combine_dfs_subject(study, "abs")
    df_np = combined_df.to_numpy()
    total_abund = df_np.sum(axis = 1)
    idxs = np.argsort(total_abund)[::-1]
    top_idxs = idxs[0 : N]
    all_taxa = study.taxa

    taxa_to_remove = [all_taxa[i].name for i in idxs[N:]]
    
    new_study = study.pop_taxa(taxa_to_remove)
    retained_taxa = [t.name for t in new_study.taxa]
    print("Retained Taxa: {}".format(retained_taxa))

    return new_study

if __name__ =="__main__":

    args = parse_arguments()
    study = md2.Study.load(args.study_file)
    top_otu_study = get_top_abundance_data(study, args.n_bugs)
    savepath = Path(args.output_loc)
    savepath.mkdir(parents=True, exist_ok=True)

    top_otu_study.save(savepath / "top_{}_otus.pkl".format(args.n_bugs))

    print("Saved study object for top {} OTUs".format(args.n_bugs))
