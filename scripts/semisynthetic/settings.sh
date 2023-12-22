# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="semisynthetic"
source ./settings.sh  # use parent settings (scripts/settings.sh)


export QPCR_NOISE_SCALE=0.01
export ALPHA_DIRICHLET=73.4451785497543  # empirical fit using Subject 2, taxa "OTU_2" at timepoint 3.0
export READ_DEPTHS=(10000 20000 30000 40000 50000)
export NUM_SAMPLE_TRIALS=5

export STUDY_PKL="${PROJECT_DIR}/datasets/gibson/preprocessed/gibson_healthy_agg_filtered.pkl"
export FIXED_MODULE_MCMC="${PROJECT_DIR}/datasets/gibson/output/mdsine2/inference/merged_studies_fixed_cluster/mcmc.pkl"
