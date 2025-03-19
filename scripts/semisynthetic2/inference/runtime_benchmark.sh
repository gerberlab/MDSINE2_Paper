#!/bin/bash
set -e
source semisynthetic2/settings.sh


function do_mdsine2_filter()
{
  source_pickle_dir=$1
  target_inference_dir=$2
  abund_threshold=$3

  mkdir -p "${target_inference_dir}"

  echo "[* runtime_benchmark.sh] Applying dataset filter."
  mdsine2 filter \
    --dataset "${source_pickle_dir}/synthetic.pkl" \
    --outfile "${target_inference_dir}/synthetic_filtered.pkl" \
    --dtype rel \
    --threshold "${abund_threshold}" \
    --min-num-consecutive 7 \
    --min-num-subjects 2 \
    --colonization-time 5

  echo "[* filter.sh] Applying dataset filter to replicates."
  python semisynthetic2/inference/mdsine2/helpers/filter_replicate.py \
    -r "${source_pickle_dir}/synthetic_replicate.pkl" \
    -l "${target_inference_dir}/synthetic_filtered.pkl" \
    -o "${target_inference_dir}/synthetic_replicate_filtered.pkl"
}

function run_inference()
{
  target_inference_dir=$1
  start_time=$(date +%s.%N)
  bash semisynthetic2/inference/mdsine2/run_mdsine2.sh "${target_inference_dir}"
  end_time=$(date +%s.%N)

  # Note: duration is stored in seconds.
  duration=$(echo "$end_time - $start_time" | bc)
  echo "${duration}" >> ${target_inference_dir}/runtime_inference.txt
}


# ======== First, run the filter.
runtime_test_dir="/data/local/semisynthetic2_runtime_test"
echo "Planning to store runtime benchmark results to ${runtime_test_dir}"

src_pickle_dir=${DATASET_DIR}/trajectory_replicate_0/perts_3/data_replicate_0/mice_4/timepoints_all/mdsine2
do_mdsine2_filter "${src_pickle_dir}" "${runtime_test_dir}/timepoints_all/abund_0_0001" 0.0001
do_mdsine2_filter "${src_pickle_dir}" "${runtime_test_dir}/timepoints_all/abund_0_0005" 0.0005
do_mdsine2_filter "${src_pickle_dir}" "${runtime_test_dir}/timepoints_all/abund_0_003" 0.003

src_pickle_dir=${DATASET_DIR}/trajectory_replicate_0/perts_3/data_replicate_0/mice_4/timepoints_thin1/mdsine2
do_mdsine2_filter "${src_pickle_dir}" "${runtime_test_dir}/timepoints_thin1/abund_0_0001" 0.0001

src_pickle_dir=${DATASET_DIR}/trajectory_replicate_0/perts_3/data_replicate_0/mice_4/timepoints_thin2/mdsine2
do_mdsine2_filter "${src_pickle_dir}" "${runtime_test_dir}/timepoints_thin2/abund_0_0001" 0.0001

src_pickle_dir=${DATASET_DIR}/trajectory_replicate_0/perts_3/data_replicate_0/mice_4/timepoints_thin3/mdsine2
do_mdsine2_filter "${src_pickle_dir}" "${runtime_test_dir}/timepoints_thin3/abund_0_0001" 0.0001

# ========= Next, run inference.
run_inference "${runtime_test_dir}/timepoints_all/abund_0_0001"
run_inference "${runtime_test_dir}/timepoints_all/abund_0_0005"
run_inference "${runtime_test_dir}/timepoints_all/abund_0_003"
run_inference "${runtime_test_dir}/timepoints_thin1/abund_0_0001"
run_inference "${runtime_test_dir}/timepoints_thin2/abund_0_0001"
run_inference "${runtime_test_dir}/timepoints_thin3/abund_0_0001"