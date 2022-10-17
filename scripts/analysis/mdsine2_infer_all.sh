#!/bin/bash
set -e
source analysis/settings.sh

all_seeds=(10 11 12 13 14 15 16 17 18 19)

for seed in ${all_seeds[@]}; do
	bash analysis/mdsine2_infer.sh $seed
done


outdir=$MDSINE2_OUT_DIR/merged_studies

echo "[*] Extracting posteriors..."
command=(mdsine2 extract-posterior -o $outdir)
for seed in ${all_seeds[@]}; do
	study_name="healthy-seed${seed}"
	mcmc=$MDSINE2_OUT_DIR/$study_name/mcmc.pkl
	command+=(-i $mcmc)
done
$command


echo "[*] Generating metadata."

metadata_file=$outdir/metadata.txt
echo "SEEDS:" > $metadata_file
for seed in ${all_seeds[@]}; do
	echo "$seed" >> $metadata_file
done
