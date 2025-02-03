#!/bin/bash
source semisynthetic2/settings.sh

traj_repl=$1
n_perts=$2
data_repl=$3
n_mice=$4
timeseries_id=$5
inference_dir=$6
require_variable "traj_repl" "${traj_repl}"
require_variable "n_perts" "${n_perts}"
require_variable "data_repl" "${data_repl}"
require_variable "n_mice" "${n_mice}"
require_variable "timeseries_id" "${timeseries_id}"
require_variable "inference_dir" "${inference_dir}"


echo "[* filter.sh] Applying dataset filter."
# filter relative abundance; relabund must exceed 0.0001 for at least 7 timepoints in 2 subjects, allowing for initial colonization time of 5 days.
if [[ "$timeseries_id" == "all" ]]; then
  bash semisynthetic2/inference/mdsine2/filter.sh ${traj_repl} ${n_perts} ${data_repl} ${n_mice} ${timeseries_id} ${inference_dir}
else
  full_timeseries_pickle_dir=${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/mice_${n_mice}/timepoints_all/mdsine2
  this_pickle_dir=${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/mice_${n_mice}/timepoints_${timeseries_id}/mdsine2

  echo "[* filter_like_densest_timeseries.sh] Applying dataset filter to replicates."
  python semisynthetic2/inference/mdsine2/helpers/filter_like_other.py \
    -i "${this_pickle_dir}/synthetic.pkl" \
    --like-other "${full_timeseries_pickle_dir}/synthetic_filtered.pkl" \
    -o "${inference_dir}/synthetic_filtered.pkl"

  echo "[* filter_like_densest_timeseries.sh] Applying dataset filter to replicates."
  python semisynthetic2/inference/mdsine2/helpers/filter_like_other.py \
    -i "${pickle_dir}/synthetic_replicate.pkl" \
    --like-other "${full_timeseries_pickle_dir}/synthetic_filtered.pkl" \
    -o "${inference_dir}/synthetic_replicate_filtered.pkl"
fi
