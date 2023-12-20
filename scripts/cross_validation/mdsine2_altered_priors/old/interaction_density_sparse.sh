#!/bin/bash
set -e
source cross_validation/settings.sh


excluded_subj=$1


seed=$((excluded_subj * 100))
inference_out_dir=$OUTPUT_DIR/mdsine2-modules-interaction-density-sparse/${excluded_subj}
mkdir -p ${inference_out_dir}

echo "[*] Running cross-validation run for MDSINE2 with sparse interaction-ind prior (Excluding subject: ${excluded_subj})"

# 870 = N * (N-1), where N=30
python cross_validation/helpers/mdsine2_loo.py \
		--input $DATASET_PKL \
		--negbin $REPLICATE_MCMC \
		--seed ${seed} \
		--burnin $MDSINE2_BURNIN \
		--n-samples $MDSINE2_SAMPLES \
		--checkpoint $MDSINE2_SAVE_EVERY \
		--basepath $inference_out_dir \
		--interaction-ind-prior "manual" \
		--interaction-ind-prior-a 870 \
		--interaction-ind-prior-b 0.5 \
		--perturbation-ind-prior $MDSINE2_PERTURBATION_INDICATOR_PRIOR \
		--log-every 100 \
		--exclude-subject $excluded_subj \
		| tee ${inference_out_dir}/log.txt
