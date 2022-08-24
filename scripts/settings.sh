# =================== DON'T MODIFY THESE (unless you really need to)
_cwd=$(pwd)
_parent="$(dirname "$_cwd")"
_this_path="${_cwd}/settings.sh"  # the location of this file.
echo "[*] Using environment settings from ${_this_path}."


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

export PROJECT_DIR=$_parent
export CLV_DIR="${PROJECT_DIR}/submodules/clv_fork"
export DATASET_DIR="${PROJECT_DIR}/datasets/${DATASET_NAME}"
export OUTPUT_DIR="${DATASET_DIR}/output"
