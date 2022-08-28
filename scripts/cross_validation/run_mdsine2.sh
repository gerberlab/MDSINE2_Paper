#!/bin/bash
set -e
source cross_validation/settings.sh


for excluded_subj in 0 1 2 3 4; do
	seed=$((excluded_subj * 100))
	inference_out_dir=$OUTPUT_DIR/mdsine2

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
			--benchmark \
			| tee ${inference_out_dir}/log.txt
