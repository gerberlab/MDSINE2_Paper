# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="semisynthetic"
source ./settings.sh  # use parent settings (scripts/settings.sh)


export QPCR_NOISE_SCALE=0.01
export ALPHA_DIRICHLET=100.0
export READ_DEPTHS=(10000 20000 30000 40000 50000)

export STUDY_PKL="${PROJECT_DIR}/datasets/gibson/preprocessed/gibson_healthy_agg_filtered.pkl"
export FIXED_MODULE_MCMC="${PROJECT_DIR}/datasets/gibson/output/mdsine2/inference/merged_studies_fixed_cluster/mcmc.pkl"
