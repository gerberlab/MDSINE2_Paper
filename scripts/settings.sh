# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export PROJECT_DIR="/mnt/d/Projects/MDSINE2_figures"
export DATASET_NAME="gibson"


# =================== DON'T MODIFY THESE (unless you really need to)
_THIS_PATH="${PROJECT_DIR}/scripts/settings.sh"  # the location of this file.
echo "[*] Using environment settings from ${_THIS_PATH}."


export DATASET_DIR="${PROJECT_DIR}/datasets/${DATASET_NAME}"
export PREPROCESS_DIR="${DATASET_DIR}/preprocessed"

# === outputs go here
export OUTPUT_DIR="${DATASET_DIR}/output"
export PHYLOGENY_OUT_DIR="${OUTPUT_DIR}/phylogeny"
export PLOTS_OUT_DIR="${OUTPUT_DIR}/plots"


require_program()
{
	command -v ${1} >/dev/null 2>&1 || {
		echo >&2 "I require ${1} but it's not installed.  Aborting.";
		exit 1;
	}
}
export require_program

