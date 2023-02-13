'''Plot the OTU abundances for each subject with the inner ASVs
'''

from pathlib import Path
import mdsine2 as md2
from mdsine2.logger import logger
import matplotlib.pyplot as plt
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--outdir', '-o', type=str, dest='outdir',
        help='This is where you want to save the parsed dataset.')
    parser.add_argument('--study', '-s', type=str, dest='study',
        help='Dataset that contains all of the information')
    parser.add_argument('--top', '-t', type=int, dest='top',
        help='Plot up to this number', default=None)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    study = md2.Study.load(args.study)

    if args.top is None:
        top = len(study.taxa)
    else:
        top = args.top

    for subj in study:
        subjpath = outdir / 'Subject_{}'.format(subj.name)
        subjpath.mkdir(exist_ok=True, parents=True)

        logger.info('Subject {}'.format(subj.name))
        for tidx, taxon in enumerate(study.taxa):
            if not isinstance(taxon, md2.OTU):
                continue
            if tidx >= top:
                break
            logger.info('taxon {}/{}'.format(taxon.idx, len(study.taxa)))

            fig = plt.figure(figsize=(10, 5))
            ax = fig.add_subplot(111)
            ax = md2.visualization.aggregate_taxa_abundances(subj=subj, agg=taxon, dtype='rel', ax=ax)
            fig.tight_layout()
            plt.savefig(subjpath / '{}.pdf'.format(taxon.name))
            plt.close()
