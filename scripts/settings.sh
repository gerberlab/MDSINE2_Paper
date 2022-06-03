# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export PROJECT_DIR="/mnt/f/MDSINE2_figures"


# =================== DON'T MODIFY THESE (unless you really need to)
_THIS_PATH="${PROJECT_DIR}/scripts/settings.sh"  # the location of this file.
echo "[*] Using environment settings from ${_THIS_PATH}."


require_program()
{
	command -v ${1} >/dev/null 2>&1 || {
		echo >&2 "I require ${1} but it's not installed.  Aborting.";
		exit 1;
	}
}

require_variable()
{
	var_name=$1
	value=$2
	if [ -z "$value" ]
	then
		echo "Environment variable \"$var_name\" is empty"
		exit 1
	fi
}

export require_program
export require_variable


require_variable "DATASET_NAME" "${DATASET_NAME}"

export DATASET_DIR="${PROJECT_DIR}/datasets/${DATASET_NAME}"
export PREPROCESS_DIR="${DATASET_DIR}/preprocessed"

export OUTPUT_DIR="${DATASET_DIR}/output"
export PHYLOGENY_OUT_DIR="${OUTPUT_DIR}/phylogeny"
export PLOTS_OUT_DIR="${OUTPUT_DIR}/plots"
