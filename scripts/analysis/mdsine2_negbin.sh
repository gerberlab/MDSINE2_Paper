#!/bin/bash

set -e
source gibson_inference/settings.sh

echo "Learning negative binomial dispersion parameters..."
echo "Output Directory: ${NEGBIN_OUT_DIR}"

mkdir -p ${NEGBIN_OUT_DIR}

mdsine2 infer-negbin \
    --input ${REPLICATE_DSET} \
    --seed 0 \
    --burnin 2000 \
    --n-samples 6000 \
    --checkpoint 200 \
    --basepath ${NEGBIN_OUT_DIR}

mdsine2 visualize-negbin \
    --chain ${REPLICATE_MCMC} \
    --output-basepath ${REPLICATE_PLOTS}

echo "Done."
