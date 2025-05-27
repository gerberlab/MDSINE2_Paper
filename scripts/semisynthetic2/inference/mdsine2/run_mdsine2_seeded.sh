#!/bin/bash
source semisynthetic2/settings.sh

inference_dir=$1
require_variable "inference_dir" "${inference_dir}"

# skip negbin qpcr fitting. (already done in run_mdsine2.sh)

# ======= Perform inference
echo "[* run_mdsine2.sh] Performing inference on dataset."
mdsine2 infer \
		--input "${inference_dir}/synthetic_filtered.pkl" \
		--negbin "${inference_dir}/synthetic-replicate/mcmc.pkl" \
		--rename-study "synthetic-seed10" \
		--seed "10" \
		--burnin "$BURNIN" \
		--n-samples "$N_SAMPLES" \
		--checkpoint "$CHECKPOINT" \
		--multiprocessing "$MULTIPROCESSING" \
		--basepath "$inference_dir" \
		--interaction-ind-prior 'weak-agnostic' \
		--perturbation-ind-prior 'weak-agnostic'
