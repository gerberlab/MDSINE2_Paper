#!/bin/bash
set -e
source cross_validation/settings.sh


excluded_subj=$1


seed=$((excluded_subj * 100))
inference_out_dir=$OUTPUT_DIR/edge_density_3/${excluded_subj}
mkdir -p ${inference_out_dir}

beta_b=50.0
echo "[*] Running cross-validation run for MDSINE2 with density (a=0.5, b=${beta_b}) (Excluding subject: ${excluded_subj})"
python cross_validation/helpers/mdsine2_loo.py \
		--input $DATASET_PKL \
		--negbin $REPLICATE_MCMC \
		--seed ${seed} \
		--burnin $MDSINE2_BURNIN \
		--n-samples $MDSINE2_SAMPLES \
		--checkpoint $MDSINE2_SAVE_EVERY \
		--basepath $inference_out_dir \
		--interaction-ind-prior "manual" \
		--interaction-ind-prior-a 0.5 \
		--interaction-ind-prior-b ${beta_b} \
		--perturbation-ind-prior $MDSINE2_PERTURBATION_INDICATOR_PRIOR \
		--log-every 100 \
		--exclude-subject $excluded_subj \
		| tee ${inference_out_dir}/log.txt
