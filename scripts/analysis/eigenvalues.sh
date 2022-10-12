#!/bin/bash
set -e
source analysis/settings.sh

require_program mdsine2
echo "Evaluating MDSINE2 learned model's eigenvalues (exclude clusters)"

# Healthy cohort
# --------------

seed=0
study="healthy-seed${seed}"
fixed_module_study="healthy-seed${seed}-fixed-cluster"

mcmc=$MDSINE2_OUT_DIR/${study}/mcmc.pkl
fixed_module_mcmc=$MDSINE2_OUT_DIR/${fixed_module_study}/mcmc.pkl

outdir=${MDSINE2_OUT_DIR}/${study}/eigenvalues
mkdir -p $outdir

# Note: index 19 = day 20
python analysis/helpers/calculate_eigenvalues.py \
		-i $MDSINE2_OUT_DIR/merged_studies \
		-o $outdir
