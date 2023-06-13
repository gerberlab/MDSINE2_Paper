#!/bin/bash
set -e
source analysis/settings.sh

echo "Running MDSINE2 model in fixed-clustering mode"
echo "Writing files to ${MDSINE2_OUT_DIR}"


study_name="merged_studies_fixed_cluster"

python analysis/helpers/mdsine2_fixed_module_multiseed.py \
		--input $HEALTHY_DSET \
		--multiseed-dir $MDSINE2_OUT_DIR/merged_studies \
		--negbin $REPLICATE_MCMC \
		--seed 0 \
		--burnin $BURNIN \
		--n-samples $N_SAMPLES \
		--checkpoint $CHECKPOINT \
		--multiprocessing $MULTIPROCESSING \
		--rename-study $study_name \
		--basepath $MDSINE2_OUT_DIR \
		--interaction-ind-prior $INTERACTION_IND_PRIOR \
		--perturbation-ind-prior $PERTURBATION_IND_PRIOR \
		| tee ${MDSINE2_OUT_DIR}/log_mdsine2_inference_fixedcluster.txt

echo "[*] Visualizing output of ${study_name}"
mdsine2 visualize-posterior \
		--chain $MDSINE2_OUT_DIR/$study_name/mcmc.pkl \
		--output-basepath $MDSINE2_OUT_DIR/$study_name/posterior

echo "Finished inference (fixed-cluster mode)."
