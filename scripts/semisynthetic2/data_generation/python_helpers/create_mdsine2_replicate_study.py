from pathlib import Path
from typing import *
import argparse

import pandas as pd
from base import generate_mouse_name
from mdsine2 import TaxaSet, Study


def load_timepoints(timepoints_text_path: Path) -> Set[int]:
    timepoints = set()
    with open(timepoints_text_path, "rt") as f:
        for line_idx, line in enumerate(f):
            line = line.strip()
            if line.startswith("#"):
                continue
            try:
                t = int(line)
            except ValueError:
                print("Unable to parse line #{} (`{}`) into number-valued timepoint. Skipping.".format(
                    line_idx + 1,
                    line
                ))
                continue

            timepoints.add(t)
    return timepoints


def main(
        metadata_path: Path,
        study_containing_taxa: Path,
        counts_path: Path,
        qpcr_path: Path,
        perts_path: Path,
        timepoints_path: Path,
        study_name: str,
        out_path: Path
):
    # First, slice the metadata by samples.
    timepoints = load_timepoints(timepoints_path)
    metadata = pd.read_csv(metadata_path, sep='\t')

    # NOTE: this line is the main difference between create_mdsine2_study.py and this script.
    subject_subset = set(f'M0-D{t}' for t in timepoints)
    metadata = metadata.loc[metadata['subject'].isin(subject_subset), :]
    if len(pd.unique(metadata['subject'])) != len(subject_subset):
        print("The following replicate subject IDs were requested but were not found in the metadata table: {}".format(
            timepoints.difference(set(pd.unique(metadata['subject'])))
        ))
        raise ValueError("The replicate subject IDs found in the dataframe did not match the requested set.")
    else:
        target_sample_ids = set(pd.unique(metadata['sampleID']))

    # Next, slice the other tables using this metadata.
    counts = pd.read_csv(counts_path, sep='\t')
    counts = counts.loc[counts["sampleID"].isin(target_sample_ids)].pivot(
        index='taxaName',
        columns='sampleID',
        values='reads'
    )

    qpcr = pd.read_csv(qpcr_path, sep='\t')
    qpcr = qpcr.loc[qpcr["sampleID"].isin(target_sample_ids)]
    try:
        perturbations = pd.read_csv(perts_path, sep='\t').set_index('name')
    except pd.errors.EmptyDataError:
        perturbations = None

    # Finally, create the pickle file.
    base_study = Study.load(str(study_containing_taxa))
    study = Study(base_study.taxa, name=study_name)

    study.parse(
        metadata=metadata.set_index('sampleID'),
        reads=counts,  # this is the result of a pivot, so index is already set properly.
        qpcr=qpcr.set_index('sampleID'),
        perturbations=perturbations,
    )
    out_path.parent.mkdir(exist_ok=True, parents=True)
    study.save(str(out_path))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", "-m", type=str, required=True)
    parser.add_argument("--taxa-study", "-t", dest="study_containing_taxa", type=str, required=True)
    parser.add_argument("--counts", "-c", type=str, required=True)
    parser.add_argument("--qpcr", "-q", type=str, required=True)
    parser.add_argument("--perturbations", "-p", type=str, required=True)
    parser.add_argument(
        "--timepoints", "-ts", type=str, required=True,
        help="The path to a text file containing a list of timepoints (one timepoint per line, as a float)"
    )
    parser.add_argument("--study-name", "-s", dest="study_name", type=str, required=True)
    parser.add_argument("--out", "-o", dest="out_path", type=str, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        metadata_path=Path(args.metadata),
        study_containing_taxa=Path(args.study_containing_taxa),
        counts_path=Path(args.counts),
        qpcr_path=Path(args.qpcr),
        perts_path=Path(args.perturbations),
        timepoints_path=Path(args.timepoints),
        study_name=args.study_name,
        out_path=Path(args.out_path),
    )
