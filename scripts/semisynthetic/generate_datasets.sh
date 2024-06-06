#!/bin/bash
set -e
source semisynthetic/settings.sh


max_read_depth=${READ_DEPTHS[-1]}
for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
  sample_seed=${trial}

  echo "[*] Sampling dataset for trial=${trial} (read depth=${max_read_depth})"
  echo "Using negbin noise model."
  raw_dataset=${DATASET_DIR}/trial_${trial}/raw/dataset.pkl
  raw_replicate=${DATASET_DIR}/trial_${trial}/raw/replicate.pkl

  python semisynthetic/generate_dataset_fwsim.py \
    --study ${STUDY_PKL} \
    --growths /data/local/MDSINE2_files/merged_studies/growth.npy \
    --interactions /data/local/MDSINE2_files/merged_studies/interactions.npy \
    --perts /data/local/MDSINE2_files/merged_studies/perturbations.npz \
    --coclust /data/local/MDSINE2_files/merged_studies_2/coclusters.npy \
    --truth-dir ${DATASET_DIR}/truth \
    --out ${raw_dataset} \
    --replicate-out ${raw_replicate} \
    --read-depth ${max_read_depth} \
    --a0 ${NEGBIN_A0} \
    --a1 ${NEGBIN_A1} \
    --qpcr-noise-scale ${QPCR_NOISE_SCALE} \
    --seed ${sample_seed} \
    --sim-dt 0.01 \
    --sim-max 1e20

  # Rarify to obtain lower read depths.
  filter_args=""
  for read_depth in "${READ_DEPTHS[@]}"; do
    rarefied_replicate=${DATASET_DIR}/trial_${trial}/read_depth_${read_depth}/replicate.pkl
    rarefied_study=${DATASET_DIR}/trial_${trial}/read_depth_${read_depth}/dataset.pkl
    filtered_study=${DATASET_DIR}/trial_${trial}/read_depth_${read_depth}/dataset_filtered.pkl
    mkdir -p ${DATASET_DIR}/trial_${trial}/read_depth_${read_depth}

    python semisynthetic/rarefy_dataset.py -s ${raw_replicate} -o ${rarefied_replicate} -r ${read_depth} --seed ${read_depth}${trial}
    python semisynthetic/rarefy_dataset.py -s ${raw_dataset} -o ${rarefied_study} -r ${read_depth} --seed ${read_depth}${trial}
    filter_args="${filter_args} -i ${rarefied_study} -o ${filtered_study}"

    # Also filter in the normal way.
    filtered_study_normal=${DATASET_DIR}/trial_${trial}/read_depth_${read_depth}/dataset_filtered_normal.pkl
    echo "Applying normal filter."
    mdsine2 filter \
        --dataset ${rarefied_study} \
        --outfile ${filtered_study_normal} \
        --dtype rel \
        --threshold 0.0001 \
        --min-num-consecutive 7 \
        --min-num-subjects 2 \
        --colonization-time 5
  done

  # filter all.
  echo "Applying jointly rarified filter."
  python semisynthetic/filter_across_rarefied.py ${filter_args} \
      --dtype rel \
      --threshold 0.0001 \
      --min-num-consecutive 7 \
      --min-num-subjects 2 \
      --colonization-time 5
done


#for read_depth in "${READ_DEPTHS[@]}"; do
#  for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
#    sample_seed=${read_depth}${trial}
#
#    echo "[*] Sampling dataset for read_depth=${read_depth} | trial=${trial}"
#    echo "Using negbin noise model."
#    python semisynthetic/generate_dataset_fwsim.py \
#      --study ${STUDY_PKL} \
#      --growths /data/local/MDSINE2_files/merged_studies/growth.npy \
#      --interactions /data/local/MDSINE2_files/merged_studies/interactions.npy \
#      --perts /data/local/MDSINE2_files/merged_studies/perturbations.npz \
#      --coclust /data/local/MDSINE2_files/merged_studies_2/coclusters.npy \
#      --cache-dir ${DATASET_DIR}/_tmp \
#      --out ${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/dataset.pkl \
#      --replicate-out ${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/replicate.pkl \
#      --read-depth ${read_depth} \
#      --a0 ${NEGBIN_A0} \
#      --a1 ${NEGBIN_A1} \
#      --qpcr-noise-scale ${QPCR_NOISE_SCALE} \
#      --seed ${sample_seed} \
#      --sim-dt 0.01 \
#      --sim-max 1e20
#
#    # ==== Healthy
#    mdsine2 filter \
#        --dataset ${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/dataset.pkl \
#        --outfile ${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/dataset_filtered.pkl \
#        --dtype rel \
#        --threshold 0.0001 \
#        --min-num-consecutive 7 \
#        --min-num-subjects 2 \
#        --colonization-time 5
#  done
#done
