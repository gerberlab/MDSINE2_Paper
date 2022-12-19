'''Reassign the taxonomy based on the consensus sequences that were 
computed with the script `preprocess.py` with RDP

Author: David Kaplan
Date: 11/30/20
MDSINE2 version: 4.0.6
'''
from pathlib import Path
import pandas as pd
import mdsine2 as md2
from mdsine2.logger import logger
import os
import argparse


def parse_rdp(fname, confidence_threshold) -> pd.DataFrame:
    '''Parse the taxonomic assignment document from RDP with a confidence
    threshold `confidence_threshold`

    Parameters
    ----------
    fname : str
        This is the name of the taxonomic assignment document from RDP
    confidence_threshold : float
        This is the minimum confidence needed for us to include it in the
        classification

    Returns
    -------
    dict( key1 -> dict ( key2 -> value ) )
        key1 : str
            OTU name
        key2 : str
            taxonomic level
        value : str
            taxonomic name
    '''
    f = open(fname, 'r')
    txt = f.read()
    f.close()

    columns = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']
    data = []
    index = []

    for i, line in enumerate(txt.split('\n')):
        if 'OTU' != line[:3]:
            continue

        splitting = line.split('%;')
        otuname = splitting[0].split(';')[0]
        splitting = splitting[1:]

        index.append(otuname)

        temp = []
        taxaed = []

        for tax_idx in range(len(splitting)):
            tax_key = columns[tax_idx]
            tax, confidence = splitting[tax_idx].replace('%', '').split(';')
            if float(confidence) > confidence_threshold:
                temp.append(tax)
                taxaed.append(tax_key)
            else:
                break

        for tax_key in columns:
            if tax_key not in taxaed:
                temp.append(md2.DEFAULT_TAXLEVEL_NAME)
        data.append(temp)

    return pd.DataFrame(data, columns=columns, index=index)


def assign_taxonomies(taxaset: md2.base.OTUTaxaSet, rdp_species_table: Path):
    logger.info(f'Parsing RDP from {rdp_species_table}')
    df = pd.read_csv(rdp_species_table, sep='\t', index_col=0)
    for otu in taxaset:
        asv = otu.components[0]
        row = df.loc[asv.name]
        otu.taxonomy = {
            'kingdom': row['Kingdom'],
            'phylum': row['Phylum'],
            'class': row['Class'],
            'order': row['Order'],
            'family': row['Family'],
            'genus': row['Genus'],
            'species': row['Species'],
            'asv': asv.name
        }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--input-study', '-i', type=str, dest='input_study_path')
    parser.add_argument('--output-study', '-o', type=str, dest='output_study_path')
    parser.add_argument('--rdp-table', '-r', type=str, dest='rdp_table',
        help='Location of RDP file')
    parser.add_argument('--confidence-threshold', '-c', type=float, dest='confidence_threshold',
        help='This is the minimum confidence required for us to use the classification')
    args = parser.parse_args()

    study = md2.Study.load(args.input_study_path)
    assign_taxonomies(study.taxa, Path(args.rdp_table))
    study.save(args.output_study_path)
