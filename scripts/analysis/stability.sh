#!/bin/bash
set -e
source analysis/settings.sh


# Healthy cohort
# --------------

module_idx_to_remove=$1
require_variable "module_idx_to_remove" $module_idx_to_remove


seed=0
study_name="healthy-seed${seed}-fixed-cluster"
mcmc=$MDSINE2_OUT_DIR/$study_name/mcmc.pkl

outdir=${MDSINE2_OUT_DIR}/${study_name}/stability
mkdir -p $outdir

python analysis/helpers/module_stability.py \
		--fixed-cluster-mcmc-path ${mcmc} \
		--study $HEALTHY_DSET \
		--module-remove-idx ${module_idx_to_remove} \
		--seed ${seed} \
		-o $outdir/stability_${module_idx_to_remove}.tsv
