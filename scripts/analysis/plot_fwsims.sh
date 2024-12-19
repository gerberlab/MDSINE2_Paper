#!/bin/bash
set -e
data_modality=$1
if [ "${data_modality}" == "healthy" ]; then
  source analysis/settings_healthy.sh
elif [ "${data_modality}" == "uc" ]; then
  source analysis/settings_uc.sh
else
  echo "data_modality argument is required and must be either 'healthy' or 'uc'. Exiting."
  exit 1
fi

require_program mdsine2
echo "Plotting forward simulations."

# Healthy cohort
# --------------

seed=0
study_name="healthy-seed${seed}"
mcmc=$MDSINE2_OUT_DIR/$study_name/mcmc.pkl

outdir=${MDSINE2_OUT_DIR}/${study_name}/trajectories
mkdir -p $outdir

for subj in 2 3 4 5; do
	subj_outdir=${outdir}/subject_${subj}
	mdsine2 forward-simulate \
			--input-mcmc ${mcmc} \
			--study ${MAIN_DSET} \
			--subject ${subj} \
			-o ${subj_outdir}/fwsims.npy \
			--plot all \
			--gibbs-subsample 100
done
