#!/bin/bash
set -e
source runtime_benchmark/settings.sh


# Runs all necessary analysis.
bash runtime_benchmark/mdsine2/infer_negbin_all.sh
bash runtime_benchmark/mdsine2/mdsine2_infer_all.sh
bash runtime_benchmark/other_methods/run_clv.sh
bash runtime_benchmark/other_methods/run_glv.sh
bash runtime_benchmark/other_methods/run_glv_ra.sh
bash runtime_benchmark/other_methods/run_lra.sh
