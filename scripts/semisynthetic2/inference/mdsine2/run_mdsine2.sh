#!/bin/bash
source semisynthetic2/settings.sh

traj_repl=$1
n_perts=$2
data_repl=$3
n_mice=$4
inference_dir=$5
require_variable "traj_repl" "${traj_repl}"
require_variable "n_perts" "${n_perts}"
require_variable "data_repl" "${data_repl}"
require_variable "n_mice" "${n_mice}"
require_variable "inference_dir" "${inference_dir}"


# ======= Fit NegBin qPCR model
echo "[*] Fitting Negative binomial model."
replicate_study=${inference_dir}/synthetic_replicate_filtered.pkl
mdsine2 infer-negbin --input "${replicate_study}" --seed "${NEGBIN_SEED}" --burnin 2000 --n-samples 6000 --checkpoint 200 --basepath "${inference_dir}"

# ======= Perform inference
echo "[*] Performing inference on dataset."
mdsine2 infer \
		--input "${inference_dir}/synthetic_filtered.pkl" \
		--negbin "${inference_dir}/synthetic-replicate/mcmc.pkl" \
		--seed "$INFERENCE_SEED" \
		--burnin "$BURNIN" \
		--n-samples "$N_SAMPLES" \
		--checkpoint "$CHECKPOINT" \
		--multiprocessing "$MULTIPROCESSING" \
		--basepath "$inference_dir" \
		--interaction-ind-prior 'weak-agnostic' \
		--perturbation-ind-prior 'weak-agnostic'
