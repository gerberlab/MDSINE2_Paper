# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="gibson"
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/functional_enrichment/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"


export PICRUST_OUTPUT_DIR="${OUTPUT_DIR}/functional_enrichment/picrust"
mkdir -p ${PICRUST_OUTPUT_DIR}
export SEQUENCES_DIR="${DATASET_DIR}/metadata_ASV/asv_sequences.fa"
export COUNTS_DIR="${DATASET_DIR}/raw_tables/counts.tsv"
export OUTPUT_REF_FILES_DIR="${OUTPUT_DIR}/functional_enrichment/ref_files"
mkdir -p "${OUTPUT_REF_FILES_DIR}"
export MDSINE2_MCMC_DIR=${OUTPUT_DIR}/mdsine2/inference/"healthy-seed0-fixed-cluster/mcmc.pkl"
