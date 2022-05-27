'''Reassign the taxonomy based on the consensus sequences that were 
computed with the script `preprocess.py` with RDP

Author: David Kaplan
Date: 11/30/20
MDSINE2 version: 4.0.6
'''
import pandas as pd
import mdsine2 as md2
from mdsine2.logger import logger
import os
import argparse


def parse_rdp(fname, confidence_threshold):
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
                temp.append(md2.pylab.base.DEFAULT_TAXLEVEL_NAME)
        data.append(temp)

    df = pd.DataFrame(data, columns=columns, index=index)
    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--rdp-table', '-r', type=str, dest='rdp_table',
        help='Location of RDP file')
    parser.add_argument('--confidence-threshold', '-c', type=float, dest='confidence_threshold',
        help='This is the minimum confidence required for us to use the classification')
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='This is where you want to save the parsed dataset.')
    args = parser.parse_args()

    logger.info('Parsing RDP')
    df = parse_rdp(fname=args.rdp_table, confidence_threshold=args.confidence_threshold)

    for dset in ['healthy', 'uc', 'replicates', 'inoculum']:
        logger.info('Replacing {}'.format(dset))
        study_fname = os.path.join(args.basepath, 'gibson_{dset}_agg.pkl'.format(dset=dset))
        study = md2.Study.load(study_fname)

        study.taxa.generate_consensus_taxonomies(df)
        study_fname = os.path.join(args.basepath, 'gibson_{dset}_agg_taxa.pkl'.format(dset=dset))
        study.save(study_fname)
