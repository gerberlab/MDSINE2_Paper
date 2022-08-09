#!/bin/bash
set -e
source runtime_benchmark/settings.sh


for (( trial = 1; trial < ${NUM_TRIALS}+1; trial++ )); do
	for n_taxa in 10 25 50 100 125 150 175 200; do
		bash runtime_benchmark/mdsine2/mdsine2_infer_nocluster.sh $n_taxa $trial
	done
done
