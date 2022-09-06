#!/bin/bash
set -e
source analysis/settings.sh

require_program mdsine2
echo "Evaluating MDSINE2 learned model's eigenvalues (exclude clusters)"

# Healthy cohort
# --------------

seed=0
study_name="healthy-seed${seed}-fixed-cluster"
mcmc=$MDSINE2_OUT_DIR/$study_name/mcmc.pkl

outdir=${MDSINE2_OUT_DIR}/${study_name}/eigenvalues
mkdir -p $outdir

# Note: index 19 = day 20
python analysis/helpers/calculate_eigenvalues.py \
		-f $mcmc \
		-o $outdir
