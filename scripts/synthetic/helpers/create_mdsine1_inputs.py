from pathlib import Path
import argparse
from typing import List

import numpy as np
import pandas as pd
import mdsine2 as md2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    # -o ${inference_out_dir} -m metadata.txt -c counts.txt -b biomass.txt
    parser.add_argument('-i', '--input_study', required=True,
                        help='<Required> The path to the input MDSINE2 study file containing the simulated data.')
    parser.add_argument('-o', '--out_dir', required=True,
                        help='<Required> The directory to save the target files to.')
    parser.add_argument('-m', '--metadata_file', required=True,
                        help='<Required> The name of the target metadata file.')
    parser.add_argument('-c', '--counts_file', required=True,
                        help='<Required> The name of the target counts file.')
    parser.add_argument('-b', '--biomass_file', required=True,
                        help='<Required> The name of the target biomass file.')
    parser.add_argument('-t', '--intervene_time_idx', required=True,
                        help='<Required> The timepoint index of the intervention.')
    return parser.parse_args()


def create_files(study: md2.Study, counts_path: Path, metadata_path: Path, biomass_path: Path, intervene_time_idxs: List[int]):
    counts_df_entries = []
    biomass_df_entries = []
    sample_ids = []

    with open(metadata_path, "wt") as metadata_file:
        print("sampleID\tisIncluded\tsubjectID\tmeasurementid\tperturbid\texptblock\tintv", file=metadata_file)
        metadata_line_idx = 1
        for subj_idx, subj in enumerate(study):
            times = subj.times
            pert_ids = np.zeros(len(times), dtype=int)
            for p_idx, pert in enumerate(study.perturbations):
                start = pert.starts[subj.name]
                end = pert.starts[subj.name]
                indices, = np.where((times >= start) & (times <= end))
                pert_ids[indices] = p_idx + 1

            for t_idx, (t, reads) in enumerate(subj.reads.items()):
                sample_id = f'SUBJ{subj_idx}_T{t}'
                counts_df_entries.append({
                    taxon.name: reads[taxa_idx]
                    for taxa_idx, taxon in enumerate(study.taxa)
                })
                sample_ids.append(sample_id)
                pert_id = pert_ids[t_idx]

                metadata_file.write(f"{sample_id}\t1\t{subj_idx}\t{t}\t{pert_id}\t1\t")
                if metadata_line_idx <= len(study.taxa):
                    metadata_file.write(str(intervene_time_idxs[metadata_line_idx - 1]))
                metadata_file.write("\n")
                metadata_line_idx += 1

                biomass_df_entries.append({
                    'mass1': subj.qpcr[t].data[0],
                    'mass2': subj.qpcr[t].data[1],
                    'mass3': subj.qpcr[t].data[2]
                })

    counts_df = pd.DataFrame(counts_df_entries)
    counts_df.index = sample_ids
    counts_df.transpose().to_csv(counts_path, index_label='#OTU ID', sep='\t')

    biomass_df = pd.DataFrame(biomass_df_entries)
    biomass_df.to_csv(biomass_path, index=False, sep='\t')


def main():
    args = parse_args()
    study = md2.Study.load(args.input_study)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    intervene_time_idxs = [args.intervene_time_idx] + [0 for _ in range(len(study.taxa) - 1)]
    create_files(study,
                 out_dir / args.counts_file,
                 out_dir / args.metadata_file,
                 out_dir / args.biomass_file,
                 intervene_time_idxs)


if __name__ == "__main__":
    main()
