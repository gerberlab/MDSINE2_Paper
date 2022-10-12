#!/bin/bash
set -e
source analysis/settings.sh


# Healthy cohort
# --------------

module_idx_to_remove=$1
require_variable "module_idx_to_remove" $module_idx_to_remove

outdir=${MDSINE2_OUT_DIR}/merged_studies
mkdir -p $outdir

sim_seed=31415

if [ "$module_idx_to_remove" == "None" ]; then
	python analysis/helpers/module_stability.py \
			-i $MDSINE2_OUT_DIR/merged_studies/stability \
			--n_module_replicates 10 \
			--study $HEALTHY_DSET \
			--seed ${sim_seed} \
			-o $outdir/stability_${module_idx_to_remove}.tsv \
			--num-trials 500
else
	python analysis/helpers/module_stability.py \
			-i $MDSINE2_OUT_DIR/merged_studies/stability \
			--n_module_replicates 10 \
			--study $HEALTHY_DSET \
			--seed ${sim_seed} \
			-o $outdir/stability_${module_idx_to_remove}.tsv \
			--num-trials 500 \
			--module-remove-idx ${module_idx_to_remove}
fi
