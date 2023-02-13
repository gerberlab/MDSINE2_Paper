# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
_cwd=$(pwd)
cd ../../..
export PROJECT_DIR="$(pwd)"
cd -

_this_path="${_cwd}/settings.sh"
echo "[*] Using settings from ${_this_path}"

# ========= Don't modify below
export LSF_DIR="${_cwd}/LSF"
export LSF_QUEUE="normal"
export LSF_CORES=4
export LSF_MEM=32000


require_program()
{
	command -v ${1} >/dev/null 2>&1 || {
		echo >&2 "I require ${1} but it's not installed.  Aborting.";
		exit 1;
	}
}
export require_program


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
export require_variable
