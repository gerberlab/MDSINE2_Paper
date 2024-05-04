import argparse
from typing import Union, List
from pathlib import Path
import numpy as np
import copy

import mdsine2 as md2
from mdsine2 import pylab as pl
from mdsine2.base import *


def main(
        study_paths: List[Path],
        out_paths: List[Path],
        abund_type: str,
        threshold: float,
        min_num_consecutive: int,
        min_num_subjects: int,
        colonization_time: int
):
    filtered_studies = consistency_filtering_multiple(
        study_array=[
            md2.Study.load(str(p))
            for p in study_paths
        ],
        abund_type=abund_type,
        threshold=threshold,
        min_num_consecutive=min_num_consecutive,
        min_num_subjects=min_num_subjects,
        colonization_time=colonization_time
    )

    for s, p1, p2 in zip(filtered_studies, study_paths, out_paths):
        s.save(str(p2))
        print("Filtered {} -> {}".format(p1, p2))


def consistency_filtering_multiple(
        study_array: List[Study],
        abund_type: str,
        threshold: Union[float, int],
        min_num_consecutive: int,
        min_num_subjects: int,
        colonization_time: Union[float, int]=None
) -> List[Study]:
    """
    An effective copy-paste (for legacy reasons, authored by @dkaplan) of md2.consistency_filtering, but modified to take multiple md2.Study objects to filter simultaneously.
    The main difference here is that the FIRST study is filtered first, then that filter is applied across all of the rest.
    Example: [study1, study2, study3] --> study1 has taxaset T = {taxa_1, ..., taxa_k}
        Then study2, study3 are also filtered using T.
    """
    if not pl.isstr(abund_type):
        raise TypeError('`abund_type` ({}) must be a str'.format(type(abund_type)))
    if abund_type not in ['raw', 'rel', 'abs']:
        raise ValueError('`abund_type` ({}) not recognized'.format(abund_type))
    if not pl.isnumeric(threshold):
        raise TypeError('`threshold` ({}) must be a numeric'.format(type(threshold)))
    if threshold <= 0:
        raise ValueError('`threshold` ({}) must be > 0'.format(threshold))
    if not pl.isint(min_num_consecutive):
        raise TypeError('`min_num_consecutive` ({}) must be an int'.format(
            type(min_num_consecutive)))
    if min_num_consecutive <= 0:
        raise ValueError('`min_num_consecutive` ({}) must be > 0'.format(min_num_consecutive))
    if colonization_time is None:
        colonization_time = 0
    if not pl.isnumeric(colonization_time):
        raise TypeError('`colonization_time` ({}) must be a numeric'.format(
            type(colonization_time)))
    if colonization_time < 0:
        raise ValueError('`colonization_time` ({}) must be >= 0'.format(colonization_time))
    if min_num_subjects is None:
        min_num_subjects = 1
    if min_num_subjects == 'all':
        min_num_subjects = len(study_array[0])
    if not pl.isint(min_num_subjects):
        raise TypeError('`min_num_subjects` ({}) must be an int'.format(
            type(min_num_subjects)))
    if min_num_subjects > len(study_array[0]) or min_num_subjects <= 0:
        raise ValueError('`min_num_subjects` ({}) value not valid'.format(min_num_subjects))

    study_copies: List[Study] = [copy.deepcopy(s) for s in study_array]

    # This bit of code determines which taxa to remove.
    first_subjset = study_copies[0]
    talley = np.zeros(len(first_subjset.taxa), dtype=int)
    for i, subj in enumerate(first_subjset):
        matrix = subj.matrix()[abund_type]
        tidx_start = None
        for tidx, t in enumerate(subj.times):
            if t >= colonization_time:
                tidx_start = tidx
                break
        if tidx_start is None:
            raise ValueError('Something went wrong')
        matrix = matrix[:, tidx_start:]

        for oidx in range(matrix.shape[0]):
            consecutive = 0
            for tidx in range(matrix.shape[1]):
                if matrix[oidx,tidx] >= threshold:
                    consecutive += 1
                else:
                    consecutive = 0
                if consecutive >= min_num_consecutive:
                    talley[oidx] += 1
                    break
    invalid_oidxs = np.where(talley < min_num_subjects)[0]
    to_delete = first_subjset.taxa.ids.order[invalid_oidxs]

    # go ahead and delte the taxa.
    for s in study_copies:
        s.pop_taxa(to_delete)

    return study_copies


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='input_paths', action='append',
                        required=True,
                        help='The dataset pickle file that you want to parse. Repeat to specify more than one.')
    parser.add_argument('--out', '-o', type=str, dest='out_paths', action='append',
                        required=True,
                        help='The desired output path. Repeat as many times as inputs were specified, and in the same order.')
    parser.add_argument('--dtype', '-d', type=str, dest='dtype',
                        required=True,
                        choices=['raw', 'rel', 'abs'],
                        help='The type of data we are using to threshold.')
    parser.add_argument('--threshold', '-t', type=float, dest='threshold',
                        required=True,
                        help='This is the threshold the taxon must pass at each timepoint')
    parser.add_argument('--min-num-consecutive', '-m', type=int, dest='min_num_consecutive',
                        required=True,
                        help='Number of consecutive timepoints to look for in a row')
    parser.add_argument('--min-num-subjects', '-s', type=int, dest='min_num_subjects',
                        required=True,
                        help='This is the minimum number of subjects this needs to be valid for.')

    parser.add_argument('--colonization-time', '-c', type=int, dest='colonization_time',
                        required=False, default=None,
                        help='This is the time we are looking after for colonization. Default to nothing')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(
        study_paths=[Path(p) for p in args.input_paths],
        out_paths=[Path(p) for p in args.out_paths],
        threshold=args.threshold,
        abund_type=args.dtype,
        colonization_time=args.colonization_time,
        min_num_subjects=args.min_num_subjects,
        min_num_consecutive=args.min_num_consecutive
    )
