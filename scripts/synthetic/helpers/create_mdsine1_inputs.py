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


def main():
    args = parse_args()
    study = md2.Study.load(args.input_study)



if __name__ == "__main__":
    main()
