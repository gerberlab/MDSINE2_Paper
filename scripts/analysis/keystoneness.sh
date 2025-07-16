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
echo "Evaluating MDSINE2 learned model's keystonenss"

# Healthy cohort
# --------------

#seed=0
#study="healthy-seed${seed}"
#fixed_module_study="healthy-seed${seed}-fixed-cluster"

#mcmc=$MDSINE2_OUT_DIR/${study}/mcmc.pkl
#fixed_module_mcmc=$MDSINE2_OUT_DIR/${fixed_module_study}/mcmc.pkl

outdir=$MDSINE2_OUT_DIR/merged_studies/keystoneness
mkdir -p $outdir

# Note: index 19 = day 20
mdsine2 evaluate-keystoneness \
		-e $MDSINE2_OUT_DIR/merged_studies \
		-s $MAIN_DSET \
		-it 19 \
		-o $outdir \
		--n-days 100 \
		--simulate-every-n 100
