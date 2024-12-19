#!/bin/bash
set -e
data_modality=$1
if [ "${data_modality}" == "healthy" ]; then
  source preprocess/settings_healthy.sh
elif [ "${data_modality}" == "uc" ]; then
  source preprocess/settings_uc.sh
else
  echo "data_modality argument is required and must be either 'healthy' or 'uc'. Exiting."
  exit 1
fi


require_program python
require_program mdsine2
echo "Performing consistency filtering over Study objects."

# Filter the OTUs using consistency filtering
mdsine2 filter \
    --dataset ${PREPROCESS_DIR}/gibson_${data_modality}_agg.pkl \
    --outfile ${PREPROCESS_DIR}/gibson_${data_modality}_agg_filtered.pkl \
    --dtype rel \
    --threshold 0.0001 \
    --min-num-consecutive 7 \
    --min-num-subjects 2 \
    --colonization-time 5

python preprocess/helpers/filter_replicates_like_other_dataset.py \
    --replicate-dataset ${PREPROCESS_DIR}/gibson_replicates_agg.pkl \
    --like-other ${PREPROCESS_DIR}/gibson_${data_modality}_agg_filtered.pkl \
    --out-path ${PREPROCESS_DIR}/gibson_replicates_agg_filtered.pkl

echo "Done."
