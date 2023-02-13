"""
Converts a .sto alignment file to .fasta.
Performs an additional task of trimming aligned columns (a column is removed if the column is only made up of gaps ---
check this behavior for the paper and make sure it makes sense!)
"""

from pathlib import Path
import argparse
from typing import List, Iterator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    return parser.parse_args()


class AlignedSeq(object):
    def __init__(self, seq_id: str, seq: str):
        self.id = seq_id
        self.seq = seq

    def append(self, seq: str):
        self.seq = ''.join([self.seq, seq])


GAP_CHARS = ['-', '~', '.']


def isgap(res):
    """Return true if given residue is a gap character
    """
    return (res in GAP_CHARS)


def read_stockholm_asvs(input_path: Path) -> Iterator[AlignedSeq]:
    objs = {}

    with open(input_path, 'rt') as in_file:
        num_asvs = 0
        for line in in_file:
            if not line.startswith("ASV_"):
                continue

            tokens = line.rstrip().split()
            asv = tokens[0]
            seq = tokens[-1].replace('.', '-')
            if asv in objs:
                objs[asv].append(seq)
            else:
                objs[asv] = AlignedSeq(asv, seq)
            num_asvs += 1
        print(f'{num_asvs} ASVs parsed.')
    yield from objs.values()


def prune_to_fasta(alns: List[AlignedSeq], out_path: Path):
    """
    Prune columns with only gaps and print alignment
    """
    # keep_cols = []
    # for i in range(len(alns[0].seq)):
    #     # deprecated: col = aln.get_column(i)
    #     col_nucs = [sr.seq[i].upper() for sr in alns]
    #     counter = Counter(col_nucs)
    #     if all([isgap(c) for c in counter.keys()]):
    #         continue
    #     keep_cols.append(i)
    #
    # print("Keeping the following columns: {}".format(', '.join([str(x + 1) for x in keep_cols])))

    with open(out_path, 'wt') as out_f:
        for aln in alns:
            print(f">{aln.id}", file=out_f)
            print(aln.seq, file=out_f)


def main():
    args = parse_args()
    asv_seqs = list(read_stockholm_asvs(Path(args.input)))
    if len({len(s.seq) for s in asv_seqs}) > 1:
        raise RuntimeError("Multiple alignment contained mismatching lengths.")
    prune_to_fasta(asv_seqs, Path(args.output))


if __name__ == "__main__":
    main()
