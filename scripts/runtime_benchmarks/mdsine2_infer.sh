#!/bin/bash
source runtime_benchmarks/settings.sh



run_inference() {
  input_pkl="${FILTER_PICKLE_DIR}/$1"
  replicate_mcmc=$2
  output_name=$3
  inference_seed=314159

  output_dir="${MDSINE2_OUT_DIR}/${output_name}"
  runtime_file=${output_dir}/${output_name}.runtime.txt
  if [ -f "$runtime_file" ]; then
    echo "runtime estimate file ${runtime_file} already exists. Skipping inference for ${input_pkl}"
    return
  else
    echo "LOGGING START/END TIMES to ${runtime_file}"
  fi

  mkdir -p "${output_dir}"
  export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/analysis/logging_to_file.ini"

  mkdir -p "${MDSINE2_OUT_DIR}"
  echo "[*] Performing MDSINE2 inference on ${input_pkl}, target output name = ${output_name}"
  echo "[*] Writing files to ${output_dir}"
  start_time=$(date +%s.%N)

  mdsine2 infer \
      --input "$input_pkl" \
      --negbin "$replicate_mcmc" \
      --seed $inference_seed \
      --burnin "$BURNIN" \
      --n-samples "$N_SAMPLES" \
      --checkpoint "$CHECKPOINT" \
      --multiprocessing "$MULTIPROCESSING" \
      --rename-study "$output_name" \
      --basepath "$MDSINE2_OUT_DIR" \
      --interaction-ind-prior "$INTERACTION_IND_PRIOR" \
      --perturbation-ind-prior "$PERTURBATION_IND_PRIOR"
  # note: mdsine2 infer auto-creates the subdirectory with the study name. (hence --basepath is set to MDSINE2_OUT_DIR)

  echo "Finished inference. Logging runtime."

  # record runtime for inference
  end_time=$(date +%s.%N)
  duration=$(echo "$end_time - $start_time" | bc)
  echo "${duration}" >> ${runtime_file}
}


#run_inference "filtered_141_n4.pkl" "${NEGBIN_OUT_DIR}/141_n4/replicates/mcmc.pkl" "141_n4"
#run_inference "filtered_141_n3.pkl" "${NEGBIN_OUT_DIR}/141_n4/replicates/mcmc.pkl" "141_n3"
#run_inference "filtered_141_n2.pkl" "${NEGBIN_OUT_DIR}/141_n4/replicates/mcmc.pkl" "141_n2"
#run_inference "filtered_141_n1.pkl" "${NEGBIN_OUT_DIR}/141_n4/replicates/mcmc.pkl" "141_n1"
run_inference "filtered_100_n4.pkl" "${NEGBIN_OUT_DIR}/100_n4/replicates/mcmc.pkl" "100_n4"
#run_inference "filtered_60_n4.pkl" "${NEGBIN_OUT_DIR}/60_n4/replicates/mcmc.pkl" "60_n4"
#run_inference "filtered_20_n4.pkl" "${NEGBIN_OUT_DIR}/20_n4/replicates/mcmc.pkl" "20_n4"
