import argparse

from Bio import SeqIO
import ete3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta_file", type=str, required=True)
    parser.add_argument("-n", "--newick_file", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    seqs = SeqIO.to_dict(SeqIO.parse(args.fasta_file, 'fasta'))
    asvs = list(seqs.keys())

    tree = ete3.Tree(args.newick_file)
    tree.prune(asvs, preserve_branch_length=True)
    tree.write(format=1, outfile=args.output)
