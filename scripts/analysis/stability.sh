#!/bin/bash
set -e
source analysis/settings.sh


# Healthy cohort
# --------------

module_idx_to_remove=$1
require_variable "module_idx_to_remove" $module_idx_to_remove


seed=0
study="healthy-seed${seed}"
fixed_module_study="healthy-seed${seed}-fixed-cluster"

mcmc=$MDSINE2_OUT_DIR/${study}/mcmc.pkl
fixed_module_mcmc=$MDSINE2_OUT_DIR/${fixed_module_study}/mcmc.pkl

outdir=${MDSINE2_OUT_DIR}/${study}/stability
mkdir -p $outdir

if [ "$module_idx_to_remove" == "None" ]; then
	python analysis/helpers/module_stability.py \
			-m $mcmc \
			-f $fixed_module_mcmc \
			--study $HEALTHY_DSET \
			--seed ${seed} \
			-o $outdir/stability_${module_idx_to_remove}.tsv
else
	python analysis/helpers/module_stability.py \
			-m $mcmc \
			-f $fixed_module_mcmc \
			--study $HEALTHY_DSET \
			--module-remove-idx ${module_idx_to_remove} \
			--seed ${seed} \
			-o $outdir/stability_${module_idx_to_remove}.tsv
fi
