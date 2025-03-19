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

echo "[*] Learning negative binomial dispersion parameters..."
echo "Output Directory: ${NEGBIN_OUT_DIR}"

mkdir -p ${NEGBIN_OUT_DIR}
mdsine2 infer-negbin \
    --input ${REPLICATE_DSET} \
    --seed 0 \
    --burnin 2000 \
    --n-samples 6000 \
    --checkpoint 200 \
    --basepath ${NEGBIN_OUT_DIR}

echo "[*] Rendering negbin visualization."
mdsine2 visualize-negbin \
    --chain ${REPLICATE_MCMC} \
    --output-basepath ${REPLICATE_PLOTS}

echo "Done."
