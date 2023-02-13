# submits the jobs for running cLV model
set -e
source settings.sh

LOG_DIR="${LSF_DIR}/logs"
LSF_DIR="${LSF_DIR}/lsf"

mkdir -p "${LSF_DIR}"
mkdir -p "${LOG_DIR}"

for seed in 10 11 12 13 14; do
    LSF_PATH="${LSF_DIR}/inference_seed_${seed}.lsf"
        # ============ Create LSF ===========
    echo
    echo "Creating lsf for inference using seed ${seed}"
    echo
    cat <<- EOFDOC > $LSF_PATH
#!/bin/bash
#BSUB -g /mdsine2_inference
#BSUB -J seed-${seed}
#BSUB -o ${LOG_DIR}/${model}_${reg}.out
#BSUB -e ${LOG_DIR}/${model}_${reg}.err
#BSUB -q ${LSF_QUEUE}

source activate mdsine2
cd ${PROJECT_DIR}/scripts

bash analysis/mdsine2_infer.sh ${seed}
EOFDOC

done
