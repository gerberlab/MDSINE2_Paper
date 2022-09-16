#!/bin/bash
set -e
source analysis/settings.sh

require_program mdsine2
echo "Evaluating MDSINE2 learned model's keystonenss"

# Healthy cohort
# --------------

seed=0
study="healthy-seed${seed}"
fixed_module_study="healthy-seed${seed}-fixed-cluster"

mcmc=$MDSINE2_OUT_DIR/${study}/mcmc.pkl
fixed_module_mcmc=$MDSINE2_OUT_DIR/${fixed_module_study}/mcmc.pkl

outdir=${MDSINE2_OUT_DIR}/${study}/keystoneness
mkdir -p $outdir

# Note: index 19 = day 20
mdsine2 evaluate-keystoneness \
		-m $mcmc \
		-f $fixed_module_mcmc \
		-s $HEALTHY_DSET \
		-it 19 \
		-o $outdir \
		--n-days 100 \
		--simulate-every-n 1
