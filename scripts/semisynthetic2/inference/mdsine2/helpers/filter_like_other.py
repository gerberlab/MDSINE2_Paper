'''Filter the `replicates` dataset so that it has identical taxa to
another dataset. We cannot pass the `replicates` dataset into a consistency
filter because there is no such thing as a consecutive timepoint in that
dataset
'''
import mdsine2 as md2
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--input-study', '-i', type=str, dest='input',
                        help='Location of the input dataset (e.g. replicate pkl)')
    parser.add_argument('--like-other', '-l', type=str, dest='other',
                        help='Location of the other dataset to filter like')
    parser.add_argument('--out-path', '-o', type=str, dest='out_path',
                        help='Location to save the output')
    args = parser.parse_args()

    input_study = md2.Study.load(args.input)
    other = md2.Study.load(args.other)

    to_delete = []
    for taxon in input_study.taxa:
        if taxon.name not in other.taxa:
            to_delete.append(taxon.name)
    input_study.pop_taxa(to_delete)
    input_study.save(args.out_path)
