#!/bin/bash

set -e
source gibson_inference/settings.sh

NEGBIN="${NEGBIN_OUT_DIR}/replicates/mcmc.pkl"
BURNIN="5000"
N_SAMPLES="15000"
CHECKPOINT="100"
MULTIPROCESSING="0"
HEALTHY_DSET="${PREPROCESS_DIR}/gibson_healthy_agg_taxa_filtered.pkl"
UC_DSET="${PREPROCESS_DIR}/gibson_uc_agg_taxa_filtered.pkl"
INTERACTION_IND_PRIOR="strong-sparse"
PERTURBATION_IND_PRIOR="weak-agnostic"

echo "Running MDSINE2 model"
echo "Writing files to ${MDSINE_OUT_DIR}"

# Healthy cohort
# --------------

for seed in 0 1; do
	study_name="healthy-seed${seed}"

	# Seed 0
	mdsine2 infer \
			--input $HEALTHY_DSET \
			--negbin $NEGBIN \
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
