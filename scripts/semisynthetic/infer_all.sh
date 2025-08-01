#!/bin/bash
set -e
source semisynthetic/settings.sh


for read_depth in "${READ_DEPTHS[@]}"; do
  for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
    bash semisynthetic/infer.sh ${read_depth} ${trial}
  done
done
