# submits the jobs for running cLV model
set -e
source settings.sh


excluded_subj=$1

require_program mdsine2
require_variable "excluded_subj" $excluded_subj

job_name="mdsine2_exclude_${excluded_subj}"
lsf_subdir="${LSF_DIR}/${job_name}"
lsf_path="${lsf_subdir}/job.lsf"
log_stdout="${lsf_subdir}/stdout.txt"
log_stderr="${lsf_subdir}/stderr.txt"
mkdir -p ${lsf_subdir}


# ============ Create LSF ===========
echo "[*] Creating lsf file ${lsf_path}"

cat <<- EOFDOC > $lsf_path
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -o ${log_stdout}
#BSUB -e ${log_stderr}
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
