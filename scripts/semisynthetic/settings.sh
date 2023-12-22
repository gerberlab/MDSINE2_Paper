# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="semisynthetic"
source ./settings.sh  # use parent settings (scripts/settings.sh)


export QPCR_NOISE_SCALE=0.01
export ALPHA_DIRICHLET=5373.2041  # empirical fit using replicate data.
export READ_DEPTHS=(600 3000 15000 75000)
export NUM_SAMPLE_TRIALS=5

export STUDY_PKL="${PROJECT_DIR}/datasets/gibson/preprocessed/gibson_healthy_agg_filtered.pkl"
export FIXED_MODULE_MCMC="${PROJECT_DIR}/datasets/gibson/output/mdsine2/inference/merged_studies_fixed_cluster/mcmc.pkl"

# inference settings
export REPLICATE_MCMC="${PROJECT_DIR}/datasets/gibson/output/mdsine2/negbin/replicates/mcmc.pkl"
export INFERENCE_SEED=31415
export BURNIN=5000
export N_SAMPLES=15000
export CHECKPOINT=100
export MULTIPROCESSING=0
export INTERACTION_IND_PRIOR="strong-sparse"
export PERTURBATION_IND_PRIOR="weak-agnostic"
