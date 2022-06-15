from pathlib import Path
import argparse

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
    return parser.parse_args()


def create_counts(study: md2.Study, counts_path: Path):
    for sample_id in study._samples:
        sid, t = study._samples[sample_id]
        if t not in study[sid].times:
            continue

        print(sid)
        print(study[sid].reads[t])


def main():
    args = parse_args()
    study = md2.Study.load(args.input_study)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)
    create_counts(study, out_dir / args.counts_file)


if __name__ == "__main__":
    main()
