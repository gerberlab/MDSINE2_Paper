import argparse
from pathlib import Path
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_fasta', type=str, required=True)
    parser.add_argument('-o', '--output_fasta', type=str, required=True)
    parser.add_argument('-s', '--start', type=int, required=True,
                        help='The first index to keep (zero-indexed).')
    parser.add_argument('-e', '--end', type=int, required=True,
                        help='The index at which to stop (zero-indexed, non-inclusive)')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    input_fasta = Path(args.input_fasta)
    output_fasta = Path(args.output_fasta)
    start_idx = args.start
    end_idx = args.end

    seqs = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
    trimmed = []
    for k, record in seqs.items():
        record.seq = record.seq[start_idx:end_idx]
        trimmed.append(record)
    SeqIO.write(trimmed, output_fasta, 'fasta')
