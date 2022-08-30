#!/bin/bash
set -e
source analysis/settings.sh

require_program mdsine2
echo "Evaluating MDSINE2 learned model's keystonenss"

# Healthy cohort
# --------------

seed=0
study_name="healthy-seed${seed}-fixed-cluster"
mcmc=$MDSINE2_OUT_DIR/$study_name/mcmc.pkl

outdir=${MDSINE2_OUT_DIR}/${study_name}/keystoneness
mkdir -p $outdir

# Note: index 19 = day 20
mdsine2 evaluate-keystoneness \
		-f $mcmc \
		-s $HEALTHY_DSET \
		-it 19 \
		-o $outdir \
		--n-days 100
