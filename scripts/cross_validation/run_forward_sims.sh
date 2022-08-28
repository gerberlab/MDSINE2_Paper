set -e
source cross_validation/settings.sh

fwsim_output_dir="${OUTPUT_DIR}/cross_validation/forward-sims-lv/"



run_fwsim()
{
	subj=$1
	model=$2
	reg=$3
	abundance_mode=$4

	echo "[*] Forward simulating, subject:${subj} model:${model}, relative abundance, elastic net"
	echo "[*] Output dir: ${fwsim_output_dir}"
  python "${CLV_DIR}/lv_forward_sims.py" \
  		-l "${CV_OUTPUT_DIR}" \
  		-r ${reg} \
  		-a ${abundance_mode} \
  		-s ${subj} \
  		-m ${model} \
  		-od "${fwsim_output_dir}" \
  		-inp "${REGRESSION_DATASET_DIR}" \
  		-is True
}


for subj in 0 1 2 3; do
	run_fwsim $subj "clv" "elastic-net" "rel"
	run_fwsim $subj "lra" "elastic-net" "rel"
	run_fwsim $subj "glv-ra" "elastic-net" "rel"
	run_fwsim $subj "glv-ra" "ridge" "rel"
	run_fwsim $subj "glv" "ridge" "abs"
	run_fwsim $subj "glv" "elastic-net" "abs"
done
