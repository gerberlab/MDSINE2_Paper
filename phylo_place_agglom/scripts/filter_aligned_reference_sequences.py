'''Filter the aligned reference sequences.

Parameters
----------
gap_length : int
    Minimum length of the gap length that we are looking at
n_rare : int
    Maximum number of sequeunces satisfying the criteria for it to be considered rare.

I used the script as:
python scripts/filter_aligned_reference_sequences.py --input-aligned tmp/rdp_combined_prefiltered_aligned.fa --input-unaligned tmp/rdp_combined_prefiltered.fa  --output tmp/rdp_combined_prefiltered_gapfiltered.fa --gap-length 1 --n-rare 5 --format fasta
    
'''

import numpy as np
import argparse

from Bio import SeqIO

def parse_args():
    '''Parse the input parameters
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-aligned', '-ia', type=str, dest='input_aligned', default=None,
        help='These are the aligned, unfiltered sequences')
    parser.add_argument('--input-unaligned', '-ir', type=str, dest='input_raw', default=None,
        help='These are the unaligned, unfiltered sequences')
    parser.add_argument('--output', '-o', type=str, dest='output_filename', default='./output.tsv',
        help='This is the output filename to save the reference names as')
    parser.add_argument('--gap-length', '-gl', type=int, dest='gap_length', default=1,
        help='How many consecutive gaps need to be there for it to be considered')
    parser.add_argument('--n-rare', '-nr', type=int, dest='n_rare', default=5,
        help='If there are less that `n-rare` nongaps at a position, then it is considered rare.')
    parser.add_argument('--format', '-f', type=str, dest='format', default='genbank',
        help='What type of format to read and write to.')
    args = parser.parse_args()
    return args

def get_sequences_to_delete(seqs, min_gap_length, max_num_seqs_in_gap, repeats_only=False):
    '''Return the sequences that have base pairs in gaps of at least
    length `min_gap_length` for `max_num_seqs_in_gap` sequences or less.

    Parameters
    ----------
    seqs : dict (str -> Record)
        maps the ID of the sequence to the record
    min_gap_length : int
        Minimum spance of the gap that we are looking at
    max_num_seqs_in_gap : int
        Maximum number of subjects that we should count as an insertion
    repeats_only : bool
        If True, only returns the sequences that come up in more than 1 distict
        insertion regions

    Returns
    -------
    list(str)
    '''
    i = 0
    for record in seqs.values():
        width = len(record.seq)
    M = np.zeros(shape=(len(seqs), width), dtype=str)

    for i, record in enumerate(seqs.values()):
        # print('here')
        M[i] = np.asarray(list(str(record.seq)))

    X = M != '.'
    Y = M != '-'
    M = X & Y

    # `M` : places where there are not gaps
    num_seqs_no_gaps = np.sum(M, axis=0)
    insertions = num_seqs_no_gaps <= max_num_seqs_in_gap

    # Get the set of ranges that satify the gap criteria
    cols = []
    curr_cols = []
    for col, val in enumerate(insertions):
        if val:
            curr_cols.append(col)
        else:
            if len(curr_cols) >= min_gap_length:
                cols.append(curr_cols)
            curr_cols = []

    # Get the sequence ids that satisfy the gap criteria
    names = list(seqs.keys())
    cnt = {}
    for curr_cols in cols:
        sums_for_rows = np.sum(M[:, curr_cols], axis=1)
        rows = np.where(sums_for_rows > 0)[0]
        curr_names = list(set([names[row] for row in rows]))

        for name in curr_names:
            if name not in cnt:
                cnt[name] = 0
            cnt[name] += 1

    names_to_ret = []
    for k,v in cnt.items():
        if repeats_only and v == 1:
            continue
        names_to_ret.append(k)

    return names_to_ret


if __name__ == '__main__':
    args = parse_args()

    seqs = SeqIO.parse(args.input_aligned, args.format)
    seqs = SeqIO.to_dict(seqs)
    seqs_to_use = SeqIO.parse(args.input_raw, args.format)
    seqs_to_use = SeqIO.to_dict(seqs_to_use)

    seqs_to_delete = set(get_sequences_to_delete(seqs, args.gap_length, args.n_rare))
    ret = []

    for k in seqs:
        try:
            v = seqs_to_use[k]
            if k not in seqs_to_delete:
                ret.append(v)
        except:
            if 'GC_RF' in k:
                # ignore
                pass
            else:
                print('{} not found'.format(k))
    
    # f = open(args.output_filename, 'w')
    # f.write('\t'.join(ret))
    # f.close()
    SeqIO.write(ret, args.output_filename, format=args.format)



