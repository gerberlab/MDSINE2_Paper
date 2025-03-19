#!/bin/bash
set -e
source runtime_benchmarks/settings.sh


run_filter() {
  threshold=$1
  target_name=$2

  echo "[* mdsine2_filter.sh] Applying dataset filter at threshold=${threshold}."
  mdsine2 filter \
    --dataset "${SOURCE_FILTERED_DSET}" \
    --dtype rel --min-num-consecutive 7 --min-num-subjects 2 --colonization-time 5 \
    --outfile "${FILTER_PICKLE_DIR}/${target_name}" \
    --threshold "${threshold}"

  echo "[* mdsine2_filter.sh] Applying dataset filter to replicates."
  python preprocess/helpers/filter_replicates_like_other_dataset.py \
    --replicate-dataset "${SOURCE_REPLICATE_DSET}" \
    --like-other "${FILTER_PICKLE_DIR}/${target_name}" \
    --out-path "${FILTER_PICKLE_DIR}/replicate_${target_name}"
}


reduce_cohort_size() {
  source_pkl=$1
  n_target_mice=$2
  target_pkl=$3

  echo "[* mdsine2_filter.sh] Reducing cohort size of ${source_pkl} to ${n_target_mice}."
  python runtime_benchmarks/helpers/reduce_cohort_size.py \
    -i "${FILTER_PICKLE_DIR}/${source_pkl}" \
    -o "${FILTER_PICKLE_DIR}/${target_pkl}" \
    -n "${n_target_mice}"
}


echo "Will write pickle files to: ${FILTER_PICKLE_DIR}"
mkdir -p ${FILTER_PICKLE_DIR}

run_filter 0.0001 "filtered_141_n4.pkl"
run_filter 0.000453 "filtered_100_n4.pkl"
run_filter 0.0016 "filtered_60_n4.pkl"
run_filter 0.023 "filtered_20_n4.pkl"

reduce_cohort_size "filtered_141_n4.pkl" 3 "filtered_141_n3.pkl"
reduce_cohort_size "filtered_141_n4.pkl" 2 "filtered_141_n2.pkl"
reduce_cohort_size "filtered_141_n4.pkl" 1 "filtered_141_n1.pkl"
