'''Make the genbank database files to make the reference tree

Parameters
----------
--bacteria-genbank : str
--archaea-genbank : str
--filtered-reference-seqs : str

I used the scipt as
python scripts/make_genbank_files_from_alignments.py --bacteria-genbank data/rdp_download_12588seqs_bacteria.gen --archaea-genbank data/rdp_download_485seqs_archaea.gen --filtered-reference-seqs tmp/rdp_combined_prefiltered_gapfiltered_aligned.fa
'''
import numpy as np
import argparse

from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bacteria-genbank', type=str, dest='bacteria_genbank',
        help='Genbank file for bacteria')
    parser.add_argument('--archaea-genbank', type=str, dest='archaea_genbank',
        help='Genbank file for archaea')
    parser.add_argument('--filtered-reference-seqs', type=str, dest='filtered_reference_seqs',
        help='Filtered reference seqs')
    args = parser.parse_args()
    return args

def main(filtered_seqs, bacteria_genbank, archaea_genbank):
    seqs_bacteria = SeqIO.parse(bacteria_genbank, format='genbank')
    # seqs_bacteria = SeqIO.to_dict(seqs_bacteria)
    seqs_archaea = SeqIO.parse(archaea_genbank, format='genbank')
    # seqs_archaea = SeqIO.to_dict(seqs_archaea)
    

    filtered_seqs = SeqIO.parse(filtered_seqs, format='fasta')
    filtered_seqs = SeqIO.to_dict(filtered_seqs)

    ret_bacteria = []
    for record in seqs_bacteria:
        if record.id in filtered_seqs:
            ret_bacteria.append(record)
    
    ret_archaea = []
    for record in seqs_archaea:
        if record.id in filtered_seqs:
            ret_archaea.append(record)

    SeqIO.write(ret_bacteria, bacteria_genbank.replace('.gen', '_filtered.gen'), format='genbank')
    SeqIO.write(ret_archaea, archaea_genbank.replace('.gen', '_filtered.gen'), format='genbank')



if __name__ == "__main__":
    
    args = parse_args()
    main(
        filtered_seqs=args.filtered_reference_seqs,
        bacteria_genbank=args.bacteria_genbank,
        archaea_genbank=args.archaea_genbank)
    