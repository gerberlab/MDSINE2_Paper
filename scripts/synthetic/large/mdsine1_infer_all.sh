#!/bin/bash
set -e
source synthetic/large/settings.sh

read_depth=25000
for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
	for noise_level in "low" "medium" "high"; do
		bash synthetic/helpers/infer_mdsine1.sh $read_depth $trial $noise_level
	done
done
