#!/bin/bash
set -e
source analysis/settings.sh

require_program mdsine2
echo "Evaluating MDSINE2 learned model's eigenvalues (exclude clusters)"

# Healthy cohort
# --------------

outdir=${MDSINE2_OUT_DIR}/merged_studies/eigenvalues
mkdir -p $outdir

# Note: index 19 = day 20
python analysis/helpers/calculate_eigenvalues.py \
		-i $MDSINE2_OUT_DIR/merged_studies \
		-o $outdir
