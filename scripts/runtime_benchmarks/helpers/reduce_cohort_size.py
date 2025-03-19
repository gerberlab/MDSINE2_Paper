import argparse
from pathlib import Path
import mdsine2 as md2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input-study', '-i', dest='input_study', required=True, type=str,
        help='The study PKL file to trim subjects from.'
    )
    parser.add_argument(
        '--out-path', '-o', dest='out_path', required=True, type=str,
        help='The path to save the resulting study object to.'
    )
    parser.add_argument(
        '--n-target-subjs', '-n', dest='n_target_subjs', required=True, type=int,
        help='The number of subjects to keep in the final study object.'
    )
    return parser.parse_args()


def main(source_study: md2.Study, out_path: Path, n_target_subjs: int):
    subj_names = [s.name for s in source_study]
    subj_names_to_delete = subj_names[n_target_subjs:]

    _ = source_study.pop_subject(subj_names_to_delete, source_study.name)
    source_study.save(str(out_path))


if __name__ == "__main__":
    args = parse_args()
    main(
        source_study=md2.Study.load(args.input_study),
        out_path=Path(args.out_path),
        n_target_subjs=args.n_target_subjs,
    )