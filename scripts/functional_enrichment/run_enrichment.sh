#!/bin/bash
set -e
source functional_enrichment/settings.sh

python functional_enrichment/helpers/run_enrichment.py --trait_abundance_table_path "${PICRUST_OUTPUT_DIR}/KO_predicted.tsv" --funit_trait_txt_file "${REF_FILES_DIR}/common/kegg_modules.txt" --mcmc_pkl_loc ${MDSINE2_MCMC_DIR} --filter1 0.75 --filter2 3 --filter3 3 --names_folder_loc "${REF_FILES_DIR}/common" --funit_type "kegg_modules" --save_loc ${PICRUST_OUTPUT_DIR}

#python functional_enrichment/helpers/run_enrichment.py --trait_abundance_table_path "${PICRUST_OUTPUT_DIR}/KO_predicted.tsv" --funit_trait_txt_file "${REF_FILES_DIR}/common/kegg_pathways.txt" --mcmc_pkl_loc ${MDSINE2_MCMC_DIR} --filter1 0.75 --filter2 1 --filter3 1 --names_folder_loc "${REF_FILES_DIR}/common" --funit_type "kegg_pathways" --save_loc ${PICRUST_OUTPUT_DIR}

python functional_enrichment/helpers/run_enrichment.py --trait_abundance_table_path "${PICRUST_OUTPUT_DIR}/EC_predicted.tsv" --funit_trait_txt_file "${REF_FILES_DIR}/common/ec_cazy.txt" --mcmc_pkl_loc ${MDSINE2_MCMC_DIR} --filter1 0.75 --filter2 3 --filter3 3 --names_folder_loc "${REF_FILES_DIR}/common" --funit_type "cazymes" --save_loc ${PICRUST_OUTPUT_DIR}
