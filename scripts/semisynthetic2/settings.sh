# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="semisynthetic2"
source ./settings.sh  # use parent settings (scripts/settings.sh)

export MAX_N_MICE=16
export QPCR_NOISE_SCALE=0.2989  #0.2989
export ALPHA_DIRICHLET=4511.6211 # old value: 5373.2041  # empirical fit using replicate data.

# from real data
export READ_DEPTH=75000
export NEGBIN_A0=4.2034E-08
export NEGBIN_A1=6.0537E-02

# Number of seeds, other data generation settings
export N_TRAJ_SEEDS=10
export N_DATA_SEEDS=3
export SIM_MAX=1e20
export SIM_DT=0.01

# Real data
export REAL_DATA_PKL="${PROJECT_DIR}/datasets/gibson/preprocessed/gibson_healthy_agg_filtered.pkl"
export FIXED_MODULE_MCMC="${PROJECT_DIR}/datasets/gibson/output/mdsine2/inference/merged_studies_fixed_cluster/mcmc.pkl"

# inference settings
export NEGBIN_SEED=314159
export INFERENCE_SEED=123456
export BURNIN=5000
export N_SAMPLES=15000
export CHECKPOINT=100
export MULTIPROCESSING=0



dataset_dir() {
  # TODO: implement the dataset hierarchy here.
  traj_seed=$0
  data_seed=$1

  # Questions to ask:
  # Variation of N_MICE: should 16 mice share randomness with 8 mice? (e.g. the first 8 mice of the 16 should always share randomness with the 8 mice set)
}
