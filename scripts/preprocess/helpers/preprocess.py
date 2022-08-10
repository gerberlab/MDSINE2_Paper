'''Preprocess (aggregate and filter) the Gibson dataset for Healthy cohort, 
Ulcerative Colitis cohort, inoculum, and replicate read datasets.

Author: David Kaplan
Date: 11/30/20
MDSINE2 version: 4.0.6

Methodology
-----------
1) Load the dataset
2) Aggregate the ASVs into OTUs using the aligned 16S v4 rRNA sequences in 
   `files/gibson_16S_rRNA_v4_seqs_aligned_filtered.fa` given a hamming-distance. 
   Once we agglomerate them together we set the sequences to the original sequence 
   (unaligned).
3) Calculate the consensus sequences
4) Rename the taxa to OTUs
5) Remove selected timepoints

The file `paper_files/preprocessing/gibson_16S_rRNA_v4_seqs_aligned_filtered.fa` 
was prepared by first aligning the ASV sequences to the reference sequeunces in the 
phylogenetic tree. Once aligned, Taxa were manually filtered out if they had poor alignment 
within the 16S rRNA v4 region. A fasta file of the Taxa removed as well as their alignments 
can be found in `paper_files/preprocessing/prefiltered_asvs.fa`. 

'''
import argparse
from Bio import SeqIO, SeqRecord, Seq
import numpy as np
from mdsine2 import Study, OTU, OTUTaxaSet

import mdsine2 as md2
from mdsine2.logger import logger
import os


def load_dataset(dataset_name: str, dataset_dir: str, max_n_species: int, sequence_file: str) -> Study:
    logger.info(f"Loaading dataset `{dataset_name}` from {dataset_dir}")
    # 1) Load the dataset
    study = md2.dataset.load_gibson(dset=dataset_name,
                                    as_df=False,
                                    species_assignment='both',
                                    load_local=dataset_dir,
                                    max_n_species=max_n_species)

    # Set the sequences for each taxon
    # Remove all taxa that are not contained in that file
    # Remove the gaps
    logger.info('Replacing sequences with the file {}'.format(sequence_file))
    seqs = SeqIO.to_dict(SeqIO.parse(sequence_file, format='fasta'))
    to_delete = []
    for taxon in study.taxa:
        if taxon.name not in seqs:
            to_delete.append(taxon.name)
    for name in to_delete:
        logger.info('Deleting {} because it was not in {}'.format(name, sequence_file))
    study.pop_taxa(to_delete)

    M = []
    for taxon in study.taxa:
        seq = list(str(seqs[taxon.name].seq))
        M.append(seq)
    M = np.asarray(M)
    n_bases = np.sum(M != '-', axis=0)
    idxs = np.where(n_bases > 0)[0]
    logger.info('There were {} valid positions out of {}.'.format(len(idxs), M.shape[1]))
    M = M[:, idxs]
    for i,taxon in enumerate(study.taxa):
        taxon.sequence = ''.join(M[i])

    return study


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='This is where you want to save the parsed dataset.')
    parser.add_argument('--hamming-distance', '-hd', type=int, dest='hamming_distance', required=True,
        help='This is the hamming radius to aggregate ASV sequences. If nothing ' \
            'is provided, then there will be no aggregation.')
    parser.add_argument('--rename-prefix', '-rp', type=str, dest='rename_prefix',
        help='This is the prefix you are renaming the aggregate taxa to. ' \
            'If nothing is provided, then they will not be renamed', default=None)
    parser.add_argument('--sequences', '-s', type=str, dest='sequences', required=True,
        help='This is the fasta file location of the aligned sequences for each taxon' \
            ' that was used for placement in the phylogenetic tree. If nothing is ' \
            'provided, then do not replace them.')
    parser.add_argument('--remove-timepoints', dest='remove_timepoints', nargs='+', default=None, 
        type=float, help='Which times to remove')
    parser.add_argument('--max-n-species', '-ms', dest='max_n_species', type=int, default=2,
        help='Maximum number of species assignments to have in the name')
    parser.add_argument('--dataset-name', '-n', dest='dataset_name', type=str, required=True,
                        help='The name of the dataset to extract.')
    parser.add_argument('--dataset-dir', '-d', dest='dataset_dir', type=str, required=True,
                        help='The directory containing the input dataset (A collection of TSV files).')

    args = parser.parse_args()
    os.makedirs(args.basepath, exist_ok=True)

    dset = args.dataset_name
    study = load_dataset(dset, args.dataset_dir, args.max_n_species, args.sequences)

    # Aggregate with specified hamming distance
    logger.info('Aggregating taxa with a hamming distance of {}'.format(args.hamming_distance))
    agg_study = md2.aggregate_items(subjset=study, hamming_dist=args.hamming_distance)

    # 3) compute consensus sequences
    if not isinstance(agg_study.taxa, OTUTaxaSet):
        raise RuntimeError("Unexpected error: Expected OTUs, but got something else. "
                           "Check internal implementation of aggregate_items().")
    agg_study.taxa.generate_consensus_seqs(threshold=0.65, noconsensus_char='N')

    # 4) Rename taxa
    if args.rename_prefix is not None:
        logger.info('Renaming taxa with prefix {}'.format(args.rename_prefix))
        agg_study.taxa.rename(prefix=args.rename_prefix, zero_based_index=False)

    # 5) Remove timepoints
    if args.remove_timepoints is not None:
        if dset in ['healthy', 'uc']:
            agg_study.pop_times(args.remove_timepoints)

    # 6) Save the study set and sequences
    logger.info("# otus: {}".format(len(agg_study.taxa)))
    agg_study.save(os.path.join(args.basepath, 'gibson_' + dset + '_agg.pkl'))
    ret = []
    for taxon in agg_study.taxa:
        ret.append(SeqRecord.SeqRecord(seq=Seq.Seq(taxon.sequence), id=taxon.name,
            description=''))
    SeqIO.write(ret, os.path.join(args.basepath, 'gibson_' + dset + '_agg.fa'), 'fasta-2line')
