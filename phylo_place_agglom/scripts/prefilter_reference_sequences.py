'''I used the script as:

python scripts/prefilter_reference_asvs.py -i data/rdp_combined.fa -o tmp/rdp_combined_prefiltered.fa --max-length 1600 --format fasta
'''
import numpy as np
import argparse

from Bio import SeqIO


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='input', default=None,
        help='These are the reference sequences that are unaligned')
    parser.add_argument('--output', '-o', type=str, dest='output_filename', default='./output.tsv',
        help='This is the output filename to save the reference names as')
    parser.add_argument('--max-length', '-n', type=int, dest='max_length', default=1600,
        help='maximum length of a sequence')
    parser.add_argument('--format', '-f', type=str, dest='format', default='fasta',
        help='What type of format to read and write to.')
    args = parser.parse_args()

    seqs = SeqIO.parse(args.input, args.format)
    seqs = SeqIO.to_dict(seqs)

    ret = []
    for k, record in seqs.items():
        if len(record.seq) <= args.max_length:
            ret.append(record)

    print('Number started:', len(seqs))
    print('Number returned:', len(ret))
    print('Number deleted:', len(seqs) - len(ret))

    
    
    SeqIO.write(ret, args.output_filename, format=args.format)