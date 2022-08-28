# submits the jobs for running cLV model
set -e
source settings.sh


excluded_subj=$1

require_program mdsine2
require_variable "excluded_subj" $excluded_subj


log_dir="${LSF_DIR}/logs"
lsf_path="${LSF_DIR}/mdsine2_${excluded_subj}.lsf"
log_stdout="${log_dir}/mdsine2_${excluded_subj}.out"
log_stderr="${log_dir}/mdsine2_${excluded_subj}.err"


mkdir -p ${LSF_DIR}
mkdir -p ${log_dir}


# ============ Create LSF ===========
echo "[*] Creating lsf file ${lsf_path}"

cat <<- EOFDOC > $lsf_path
#!/bin/bash
#BSUB -J mdsine2_${excluded_subj}
#BSUB -o $log_stdout
#BSUB -e $log_stderr
#BSUB -q ${LSF_QUEUE}

echo '---PROCESS RESOURCE LIMITS---'
ulimit -a
echo '---LSF Parameters:---'
printenv | grep '^LSF'
echo '---LSB Parameters:---'
printenv | grep '^LSB'

cd ${PROJECT_DIR}/scripts
bash cross_validation/run_mdsine2.sh ${excluded_subj}
EOFDOC

done
