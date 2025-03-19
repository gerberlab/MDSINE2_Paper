#!/bin/bash
source runtime_benchmarks/settings.sh


fit_replicate_data() {
  replicate_pkl="${FILTER_PICKLE_DIR}/$1"
  outdir="${NEGBIN_OUT_DIR}/$2"

  mdsine2 infer-negbin \
    --input ${replicate_pkl} \
    --seed 0 \
    --burnin 2000 \
    --n-samples 6000 \
    --checkpoint 200 \
    --basepath ${outdir}
}

echo "[*] Learning negative binomial dispersion parameters..."
echo "Calibration Output Directory: ${NEGBIN_OUT_DIR}"
mkdir -p "${NEGBIN_OUT_DIR}"


fit_replicate_data "replicate_filtered_141_n4.pkl" "141_n4"
fit_replicate_data "replicate_filtered_100_n4.pkl" "100_n4"
fit_replicate_data "replicate_filtered_60_n4.pkl" "60_n4"
fit_replicate_data "replicate_filtered_20_n4.pkl" "20_n4"
