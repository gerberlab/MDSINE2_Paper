#!/bin/bash
set -e
source cross_validation/settings.sh


excluded_subj=$1

N_THREADS=16  # alienware has 16 cores=32 threads, we want to use 8 cores=16 threads.
export MKL_NUM_THREADS=${N_THREADS}
export NUMEXPR_NUM_THREADS=${N_THREADS}
export OMP_NUM_THREADS=${N_THREADS}
echo "USING ${N_THREADS} THREADS"


seed=$((excluded_subj * 100))
inference_out_dir=$OUTPUT_DIR/mdsine2-noqpcr-modules/${excluded_subj}
mkdir -p ${inference_out_dir}
echo "OUTPUT DIR: ${inference_out_dir}"

echo "[*] Running cross-validation run for MDSINE2 (Excluding subject: ${excluded_subj})"
python cross_validation/helpers/mdsine2_loo.py \
		--input $DATASET_PKL \
		--negbin $REPLICATE_MCMC \
		--seed ${seed} \
		--burnin $MDSINE2_BURNIN \
		--n-samples $MDSINE2_SAMPLES \
		--checkpoint $MDSINE2_SAVE_EVERY \
		--basepath $inference_out_dir \
		--interaction-ind-prior $MDSINE2_INTERACTION_INDICATOR_PRIOR \
		--perturbation-ind-prior $MDSINE2_PERTURBATION_INDICATOR_PRIOR \
		--log-every 100 \
		--exclude-subject $excluded_subj \
		--no-qpcr \
		| tee ${inference_out_dir}/log.txt
