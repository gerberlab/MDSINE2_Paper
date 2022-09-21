# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="gibson"
set -e 

source preprocess/settings.sh

export PREPROCESSED_ALL_DIR="${DATASET_DIR}/preprocessed_all" 
export FIGURES_OUTPUT_DIR="${OUTPUT_DIR}/figures"
export ENRICHMENT_OUTPUT_DIR="${OUTPUT_DIR}/functional_enrichment"


