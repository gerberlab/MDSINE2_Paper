#!/bin/bash
set -e
source synthetic/small/settings.sh


require_program python

read_depth=25000
for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
	seed=$trial
	python synthetic/helpers/create_datasets.py \
	-i ${GLV_PARAMS} \
	-t ${TIME_POINTS} \
	-n ${COHORT_SIZE} \
	-o ${DATASET_DIR}/data/trial_${trial} \
	-s ${seed} \
	--num_qpcr_triplicates ${QPCR_TRIPLICATES} \
	--process_var ${PROCESS_VAR} \
	-dt ${SIMULATION_DT} \
	--read_depth $read_depth \
	--low_noise ${LOW_NOISE_SCALE} \
	--medium_noise ${MEDIUM_NOISE_SCALE} \
	--high_noise ${HIGH_NOISE_SCALE} \
	--intervene_day 10.0
done
