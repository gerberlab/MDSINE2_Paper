#!/bin/bash
set -e
data_modality=$1
if [ "${data_modality}" == "healthy" ]; then
  source analysis/settings_healthy.sh
  all_seeds=(10 11 12 13 14 15 16 17 19 20)   # seeds for healthy dataset (others attempted encounter a numeric precision crash)
elif [ "${data_modality}" == "uc" ]; then
  source analysis/settings_uc.sh
  all_seeds=(10 12 13 14 15 17 19 26 30 31)   # seeds for uc dataset (others attempted encounter a numeric precision crash)
else
  echo "data_modality argument is required and must be either 'healthy' or 'uc'. Exiting."
  exit 1
fi


## Perform inference.
#for seed in ${all_seeds[@]}; do
#	bash analysis/mdsine2_infer.sh ${data_modality} $seed
#done



outdir=$MDSINE2_OUT_DIR/merged_studies
# TMP CODE: copy from dan's workflow
#echo " copying dan's workflow."
#for seed in ${all_seeds[@]}; do
#	study_name="${data_modality}-seed${seed}"
#  rsync -av /data/bwh-comppath-full/gibsonlab/dan/MDSINE2_UC/MDSINE2_Paper/datasets/gibson/output/mdsine2/inference/${study_name} $MDSINE2_OUT_DIR
#done

echo "[*] Extracting posteriors..."
# Extract posterior values from pickle.
# Build the command in a for-loop (using multiple -i argument instances).
command="mdsine2 extract-posterior -o ${outdir}"
for seed in ${all_seeds[@]}; do
	study_name="${data_modality}-seed${seed}"
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
