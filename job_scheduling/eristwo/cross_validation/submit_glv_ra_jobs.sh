# submits the jobs for running cLV model
set -e 
source settings.sh 

abund_type="rel"
model="glv-ra"
LOG_DIR="${LSF_DIR}/logs_${model}"
LSF_DIR="${LSF_DIR}/lsf_${model}"

mkdir -p "${LSF_DIR}"
mkdir -p "${LOG_DIR}"

regression_types=("ridge" "elastic-net")

for reg in "${regression_types[@]}"; do
    LSF_PATH="${LSF_DIR}/${model}_${reg}_${abund_type}.lsf"
        # ============ Create LSF ===========
    echo 
    echo "Creating lsf for ${model}, ${reg} regression and ${abund_type} abundance"
    echo
    cat <<- EOFDOC > $LSF_PATH
#!/bin/bash
#BSUB -J ${model}_${reg}
#BSUB -o ${LOG_DIR}/${model}_${reg}.out
#BSUB -e ${LOG_DIR}/${model}_${reg}.err
#BSUB -q ${LSF_QUEUE}

source activate mdsine2
cd ${PROJECT_DIR}/scripts

bash cv_comparator_methods/run_${model}.sh $reg
EOFDOC

done
