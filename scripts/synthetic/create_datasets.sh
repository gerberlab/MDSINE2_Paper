#!/bin/bash
set -e
source synthetic/settings.sh


require_program python
cd synthetic/helpers


for read_depth in 1000 25000; do
	for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
		seed=$trial
		python create_datasets.py \
		-i ${GLV_PARAMS} \
		-t ${TIME_POINTS} \
		-n ${COHORT_SIZE} \
		-o ${DATASET_DIR}/data/reads_${read_depth}/trial_${trial} \
		-s ${seed} \
		--process_var ${PROCESS_VAR} \
		-dt ${SIMULATION_DT} \
		--read_depth $read_depth \
		--low_noise ${LOW_NOISE_SCALE} \
		--medium_noise ${MEDIUM_NOISE_SCALE} \
		--high_noise ${HIGH_NOISE_SCALE}
	done
done
