# submits the jobs for running cLV model
set -e
source settings.sh


model=$1
regression_type=$2


require_program mdsine2
require_variable "model" $model
require_variable "regression_type" $regression_type


log_dir="${LSF_DIR}/logs"
lsf_path="${LSF_DIR}/${model}_${regression_type}.lsf"
log_stdout="${log_dir}/${model}_${regression_type}.out"
log_stderr="${log_dir}/${model}_${regression_type}.err"


mkdir -p ${LSF_DIR}
mkdir -p ${log_dir}


# ============ Create LSF ===========
echo "[*] Creating lsf file ${lsf_path}"

cat <<- EOFDOC > $lsf_path
#!/bin/bash
#BSUB -J ${model}_${regression_type}
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
bash cross_validation/run_regression.sh ${model} ${regression_type}
EOFDOC

done
