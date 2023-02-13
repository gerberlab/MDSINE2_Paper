import argparse
from Bio import Phylo


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    print('\n\nMaking newick tree')
    Phylo.convert(args.input, 'phyloxml', args.output, 'newick')
