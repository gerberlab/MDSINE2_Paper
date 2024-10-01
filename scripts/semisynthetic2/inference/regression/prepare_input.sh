#!/bin/bash
set -e
source semisynthetic2/settings.sh


out_dir=$1
mdsine2_pkl=$2
require_variable "out_dir" "${out_dir}"
require_variable "mdsine2_pkl" "${mdsine2_pkl}"


set -e
echo "[* prepare_input.sh] Creating CLV files in ${out_dir}"
python semisynthetic/regression/create_clv_inputs.py -i "${mdsine2_pkl}" -o "${out_dir}"
set +e
