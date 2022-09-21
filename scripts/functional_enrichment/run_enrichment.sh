#!/bin/bash
set -e
source functional_enrichment/settings.sh

python functional_enrichment/helpers/run_enrichment.py --trait_abundance_table_path "${PICRUST_OUTPUT_DIR}/KO_predicted.tsv" --funit_trait_txt_file "${OUTPUT_REF_FILES_DIR}/kegg_modules_ko.txt" --mcmc_pkl_loc ${MDSINE2_MCMC_DIR} --filter1 0.75 --filter2 3 --filter3 3 --names_file_loc "${OUTPUT_REF_FILES_DIR}/modules_names.csv" --funit_type "kegg_modules" --save_loc ${PICRUST_OUTPUT_DIR}  --pathway_module_pkl_file "${OUTPUT_REF_FILES_DIR}/module_details.pkl"

python functional_enrichment/helpers/run_enrichment.py --trait_abundance_table_path "${PICRUST_OUTPUT_DIR}/KO_predicted.tsv" --funit_trait_txt_file "${OUTPUT_REF_FILES_DIR}/kegg_pathways_ko.txt" --mcmc_pkl_loc ${MDSINE2_MCMC_DIR} --filter1 0.75 --filter2 3 --filter3 3 --names_file_loc "${OUTPUT_REF_FILES_DIR}/pathways_names.csv" --funit_type "kegg_pathways" --save_loc ${PICRUST_OUTPUT_DIR}  --pathway_module_pkl_file "${OUTPUT_REF_FILES_DIR}/pathway_details.pkl"

python functional_enrichment/helpers/run_enrichment.py --trait_abundance_table_path "${PICRUST_OUTPUT_DIR}/EC_predicted.tsv" --funit_trait_txt_file "${OUTPUT_REF_FILES_DIR}/ec_cazy.txt" --mcmc_pkl_loc ${MDSINE2_MCMC_DIR} --filter1 0.75 --filter2 3 --filter3 3 --names_file_loc "${OUTPUT_REF_FILES_DIR}/cazyme_names.csv" --funit_type "cazymes" --save_loc ${PICRUST_OUTPUT_DIR}  
