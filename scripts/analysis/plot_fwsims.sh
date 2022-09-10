#!/bin/bash
set -e
source analysis/settings.sh

require_program mdsine2
echo "Evaluating MDSINE2 learned model's eigenvalues (exclude clusters)"

# Healthy cohort
# --------------

seed=0
study_name="healthy-seed${seed}"
mcmc=$MDSINE2_OUT_DIR/$study_name/mcmc.pkl

outdir=${MDSINE2_OUT_DIR}/${study_name}/trajectories
mkdir -p $outdir

for subj in 2 3 4 5; do
	mdsine2 forward-simulate \
			--input-mcmc ${mcmc} \
			--study ${HEALTHY_DSET} \
			--subject ${subj} \
			-o ${outdir}/fwsims.npy \
			--plot all
done
