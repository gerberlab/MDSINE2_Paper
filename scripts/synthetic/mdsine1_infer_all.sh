#!/bin/bash
set -e
source synthetic/settings.sh

for read_depth in 25000; do
	for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
		for noise_level in "low" "medium" "high"; do
			bash synthetic/helpers/infer_mdsine1.sh $read_depth $trial $noise_level
		done
	done
done
