import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_fasta', type=str, required=True,
                        help='The path to the aligned Fasta file.')
    parser.add_argument('-t', '--threshold', type=int, required=True,
                        help='A length threshold for non-gapped nucleotides. Any ASV whose post-aligned sequence '
                             'contains fewer than this many nucleotides will be excluded.')
    parser.add_argument('-o', '--output_fasta', type=str, required=True,
                        help='The target output path to write to (in fasta format).')
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.input_fasta, 'r') as in_file, open(args.output_fasta, 'w') as out_file:
        for record in SeqIO.parse(in_file, "fasta"):
            record.seq = Seq(str(record.seq).replace("-", ""))
            if len(record.seq) < args.threshold:
                print(f"Truncating record {record.id} (length = {len(record.seq)})")
            else:
                SeqIO.write(record, out_file, "fasta")


if __name__ == "__main__":
    main()
