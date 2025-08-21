#!/bin/bash
set -e
source cross_validation/settings.sh


excluded_subj=$1


seed=$((excluded_subj * 100))
inference_out_dir=$OUTPUT_DIR/other_priors/pert_var_4/${excluded_subj}
mkdir -p ${inference_out_dir}

pert_var_mean=1e6  # skip 1e4; that's the default setting.
echo "[*] Running cross-validation run for MDSINE2 with perturbation variance mean=${pert_var_mean} (Excluding subject: ${excluded_subj})"
python cross_validation/helpers/mdsine2_loo_sics.py \
		--input $DATASET_PKL \
		--negbin $REPLICATE_MCMC \
		--seed ${seed} \
		--burnin $MDSINE2_BURNIN \
		--n-samples $MDSINE2_SAMPLES \
		--checkpoint $MDSINE2_SAVE_EVERY \
		--basepath $inference_out_dir \
		--interaction-ind-prior $MDSINE2_INTERACTION_INDICATOR_PRIOR \
		--perturbation-ind-prior $MDSINE2_PERTURBATION_INDICATOR_PRIOR \
		--pert-var-mean ${pert_var_mean} \
		--log-every 100 \
		--exclude-subject $excluded_subj \
		| tee ${inference_out_dir}/log.txt
