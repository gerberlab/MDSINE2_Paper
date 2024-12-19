#!/bin/bash
set -e
data_modality=$1
if [ "${data_modality}" == "healthy" ]; then
  source analysis/settings_healthy.sh
elif [ "${data_modality}" == "uc" ]; then
  source analysis/settings_uc.sh
else
  echo "data_modality argument is required and must be either 'healthy' or 'uc'. Exiting."
  exit 1
fi

all_seeds=(10 11 12 13 14 15 16 17 19 20)

# Perform inference.
for seed in ${all_seeds[@]}; do
	bash analysis/mdsine2_infer.sh ${data_modality} $seed
done


outdir=$MDSINE2_OUT_DIR/merged_studies

echo "[*] Extracting posteriors..."
# Extract posterior values from pickle.
# Build the command in a for-loop (using multiple -i argument instances).
command="mdsine2 extract-posterior -o ${outdir}"
for seed in ${all_seeds[@]}; do
	study_name="healthy-seed${seed}"
	mcmc=$MDSINE2_OUT_DIR/$study_name/mcmc.pkl
	command="${command} -i ${mcmc}"
done
$command


echo "[*] Generating metadata."
# Generate metadata file (which seeds were run?)
metadata_file=$outdir/metadata.txt
echo "SEEDS:" > $metadata_file
for seed in ${all_seeds[@]}; do
	echo "$seed" >> $metadata_file
done
