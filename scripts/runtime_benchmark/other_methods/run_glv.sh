#!/bin/bash
set -e
source runtime_benchmark/settings.sh

regression_types=("elastic-net" "ridge")

for n_taxa in 10 25 50 100; do
	for reg in "${regression_types[@]}"; do
		bash runtime_benchmark/other_methods/run_single_method.sh ${n_taxa} "glv" ${reg}
	done 
done 

