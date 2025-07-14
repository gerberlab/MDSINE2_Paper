_cwd=$(pwd)
_parent="$(dirname "$_cwd")"
PROJECT_DIR=$_parent
echo $PROJECT_DIR
DATASET_DIR="${PROJECT_DIR}/datasets/gibson"

mkdir -p gibson_full_timepoints


for data_modality in "healthy" "uc"; do
  echo "[*] Working in modality=${data_modality}"
  asv_aln_fasta=${DATASET_DIR}/${data_modality}/metadata_ASV/aligned_asvs.fa
  prefilt_aln_fasta=${DATASET_DIR}/${data_modality}/metadata_ASV/prefiltered_asvs.fa
  mkdir -p ./gibson_full_timepoints/${data_modality}

  for dataset in "${data_modality}" replicates inoculum; do
    echo "[**] Extracting dataset: ${dataset}"
    python ${PROJECT_DIR}/scripts/preprocess/helpers/preprocess.py \
        --hamming-distance 0 \
        --sequences $prefilt_aln_fasta \
        --max-n-species 2 \
        --dataset-name ${dataset} \
        --dataset-dir ${DATASET_DIR}/${data_modality}/raw_tables \
        --output-prefix "gibson_full_timepoints/${data_modality}/gibson_${dataset}_agg_full_timepoints" \
        --trim-option "ALL_GAPS" \
        --sort-order "MIN_ASV_IDX" \
        --naming-scheme "MIN_ASV_LABEL"
  done
done
