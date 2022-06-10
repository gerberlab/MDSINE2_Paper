#!/bin/bash
set -e
source runtime_benchmark/settings.sh


for n_taxa in 10 25 50 100; do
	for (( trial = 1; trial < ${NUM_TRIALS}+1; trial++ )); do
		bash runtime_benchmark/mdsine2/mdsine2_infer.sh $n_taxa $trial
	done
done
