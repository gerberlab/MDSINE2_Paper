set -e
source cv_comparator_methods/settings.sh

subjects=(0 1 2 3)
model_elastic_rel=("clv" "lra" "glv-ra")
model_ridge_rel=("glv-ra")

#for elastic net regression relative abundance
for subj in "${subjects[@]}"; do
    for model in "${model_elastic_rel[@]}"; do
        echo "Forward simulating, subject:${subj} model:${model}, relative abundance, elastic net"
        python "${CLV_DIR}/lv_forward_sims.py" -l "${CV_OUTPUT_DIR}" -r "elastic-net" -a "rel" -s ${subj} -m ${model} -od "${FW_SIM_OUTPUT_DIR}" -inp "${INPUT_DATASET_DIR}" -is True
        echo
    done
done

#for ridge regression relative abundance
for subj in "${subjects[@]}"; do
    for model in "${model_ridge_rel[@]}"; do
        echo "Forward simulating, subject:${subj} model:${model}, relative abundance, ridge"
        python "${CLV_DIR}/lv_forward_sims.py" -l "${CV_OUTPUT_DIR}" -r "ridge" -a "rel" -s ${subj} -m ${model} -od "${FW_SIM_OUTPUT_DIR}" -inp "${INPUT_DATASET_DIR}" -is True
    done
done

regression_types=("ridge" "elastic-net")
for subj in "${subjects[@]}"; do
    for reg in "${regression_types[@]}"; do
        echo "Forward simulating, subject:${subj} regression:${reg}, absolute abundance, glv"
        python "${CLV_DIR}/lv_forward_sims.py" -l "${CV_OUTPUT_DIR}" -r ${reg} -a "abs" -s ${subj} -m "glv" -od "${FW_SIM_OUTPUT_DIR}" -inp "${INPUT_DATASET_DIR}" -is True
    done
done

