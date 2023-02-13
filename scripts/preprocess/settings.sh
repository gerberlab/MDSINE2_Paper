# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="gibson"
source ./settings.sh  # use parent settings (scripts/settings.sh)
export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/preprocess/logging.ini"

export REFERENCE_RDP_DIR="${DATASET_DIR}/metadata_REFERENCE"
export READS_DIR="${DATASET_DIR}/reads"
export ASV_FASTA="${DATASET_DIR}/metadata_ASV/asv_sequences.fa"
export OTU_FASTA="${DATASET_DIR}/metadata_OTU/otu_sequences.fa"
export PREPROCESS_DIR="${DATASET_DIR}/preprocessed"

export PREFILT_LEN=250
