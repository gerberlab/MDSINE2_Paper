#!/bin/bash
set -e
source analysis/settings.sh

echo "Running MDSINE2 model"
echo "Writing files to ${MDSINE2_OUT_DIR}"

# Healthy cohort
# --------------

for seed in 0 1; do
	study_name="healthy-seed${seed}"

	# Seed 0
	echo "[*] Performing MDSINE2 inference on ${study_name}"
	mdsine2 infer \
			--input $HEALTHY_DSET \
			--negbin $REPLICATE_MCMC \
			--seed $seed \
			--burnin $BURNIN \
			--n-samples $N_SAMPLES \
			--checkpoint $CHECKPOINT \
			--multiprocessing $MULTIPROCESSING \
			--rename-study $study_name \
			--basepath $MDSINE2_OUT_DIR \
			--interaction-ind-prior $INTERACTION_IND_PRIOR \
			--perturbation-ind-prior $PERTURBATION_IND_PRIOR

	mdsine2 visualize-posterior \
			--chain $MDSINE2_OUT_DIR/$study_name/mcmc.pkl \
			--output-basepath $MDSINE2_OUT_DIR/$study_name/posterior

	echo "Finished Healthy (seed ${seed})."
done
