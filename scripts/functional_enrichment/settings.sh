# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="gibson"
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/functional_enrichment/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"


export PICRUST_OUTPUT_DIR="${OUTPUT_DIR}/functional_enrichment"
mkdir -p ${PICRUST_OUTPUT_DIR}
export SEQUENCES_DIR="${DATASET_DIR}/metadata_ASV/asv_sequences.fa"
export COUNTS_DIR="${DATASET_DIR}/raw_tables/counts.tsv"
export REF_FILES_DIR="${DATASET_DIR}/functional_enrichment"
export MDSINE2_MCMC_DIR=${OUTPUT_DIR}/mdsine2/inference/"healthy-seed0-fixed-cluster/mcmc.pkl"
