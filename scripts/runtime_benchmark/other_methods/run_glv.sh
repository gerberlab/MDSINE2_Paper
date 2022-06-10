#!/bin/bash
set -e
source runtime_benchmark/settings.sh


for n_taxa in 10 25 50 100; do
	for (( trial = 1; trial < ${NUM_TRIALS}+1; trial++ )); do
		bash runtime_benchmark/other_methods/run_single_method.sh ${n_taxa} "glv" "elastic-net"
		bash runtime_benchmark/other_methods/run_single_method.sh ${n_taxa} "glv" "ridge"
	done
done
