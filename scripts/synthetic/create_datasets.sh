#!/bin/bash
set -e
source synthetic/settings.sh


require_program python
cd synthetic/helpers


for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
	python create_datasets.py \
	-i ${GLV_PARAMS} \
	-t ${TIME_POINTS} \
	-n ${COHORT_SIZE} \
	-o ${DATASET_DIR}/datasets/trial_${seed} \
	-s ${seed} \
	--process_var ${PROCESS_VAR} \
	-dt ${SIMULATION_DT} \
	-a0 ${NEGBIN_A0} \
	-a1 ${NEGBIN_A1} \
	--read_depth ${READ_DEPTH} \
	--low_noise ${LOW_NOISE_SCALE} \
	--medium_noise ${MEDIUM_NOISE_SCALE} \
	--high_noise ${HIGH_NOISE_SCALE}
done
