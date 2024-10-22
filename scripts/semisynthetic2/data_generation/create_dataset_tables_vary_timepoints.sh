#!/bin/bash
set -e
source semisynthetic2/settings.sh


n_mice=4  # matches real data
n_perts=3  # matches real data

#declare -a timeseries_ids=("all" "thin1" "thin2" "thin3" "A" "B" "C")
declare -a timeseries_ids=("A" "B" "C" "all" "thin1" "thin2" "thin3")


for timeseries_id in "${timeseries_ids[@]}"; do
  for (( traj_repl = 0; traj_repl < ${N_TRAJ_SEEDS}; ++traj_repl)); do
    for (( data_repl = 0; data_repl < ${N_DATA_SEEDS}; ++data_repl )); do
      pickle_outdir=${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/mice_${n_mice}/timepoints_${timeseries_id}/mdsine2
      mkdir -p "${pickle_outdir}"

      dataset_timepoints_file="semisynthetic2/data_generation/files/inference_timepoints_${timeseries_id}.txt"
      replicate_timepoints_file="semisynthetic2/data_generation/files/replicate_timepoints.txt"

      # dataset for inference
      echo "[* create_dataset_tables_vary_mice.sh] Creating inference pickle file for [traj_repl=${traj_repl} | n_perts=${n_perts} | data_repl=${data_repl} | n_mice=${n_mice} | timeseries=${timeseries_id}]"
      python semisynthetic2/data_generation/python_helpers/create_mdsine2_study.py \
      -m "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/tables_raw/metadata.tsv" \
      -t "${REAL_DATA_PKL}" \
      -c "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/tables_raw/counts.tsv" \
      -q "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/tables_raw/qpcr.tsv" \
      -p "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/trajectories/perturbations.tsv" \
      -n "$n_mice" \
      -ts "$dataset_timepoints_file" \
      -s "synthetic" \
      -o "${pickle_outdir}/synthetic.pkl"

      # replicate dataset for calibration
      echo "[* create_dataset_tables_vary_mice.sh] Creating replicate pickle file for [traj_repl=${traj_repl} | n_perts=${n_perts} | data_repl=${data_repl} | n_mice=${n_mice} | timeseries=${timeseries_id}]"
      python semisynthetic2/data_generation/python_helpers/create_mdsine2_replicate_study.py \
      -m "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/tables_raw/replicate_metadata.tsv" \
      -t "${REAL_DATA_PKL}" \
      -c "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/tables_raw/replicate_counts.tsv" \
      -q "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/tables_raw/replicate_qpcr.tsv" \
      -p "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/trajectories/perturbations.tsv" \
      -ts "$replicate_timepoints_file" \
      -s "synthetic-replicate" \
      -o "${pickle_outdir}/synthetic_replicate.pkl"
    done
  done
done