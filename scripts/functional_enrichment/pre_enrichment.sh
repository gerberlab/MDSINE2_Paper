#!/bin/sh
set -e
source functional_enrichment/settings.sh

echo "Scraping information on Kegg Orthologs. Saving output to ${OUTPUT_REF_FILES_DIR}"
python functional_enrichment/helpers/KO_scraper.py --save_loc "${OUTPUT_REF_FILES_DIR}" --output_name "ko_details"

echo

echo "Scraping information on Kegg Modules. Saving output to ${OUTPUT_REF_FILES_DIR}"

echo

python functional_enrichment/helpers/pathway_module_scraper.py --ko_details_dir "${OUTPUT_REF_FILES_DIR}/ko_details.pkl" --save_loc ${OUTPUT_REF_FILES_DIR} --kegg_type "module"

echo

echo "Scraping information on Kegg Pathways. Saving output to ${OUTPUT_REF_FILES_DIR}"

python functional_enrichment/helpers/pathway_module_scraper.py --ko_details_dir "${OUTPUT_REF_FILES_DIR}/ko_details.pkl" --save_loc ${OUTPUT_REF_FILES_DIR} --kegg_type "pathway"

#python functional_enrichment/helpers/obtain_ko_category.py --funit_ko_file "${REF_FILES_DIR}/common/kegg_modules.txt" --save_loc "${REF_FILES_DIR}/common" --save_name "ko_category.txt"

#python functional_enrichment/helpers/obtain_kegg_category.py --funit_transit_txt_file_path "${REF_FILES_DIR}/common/kegg_modules.txt" --save_loc "${REF_FILES_DIR}/common" --ko_category_file "${REF_FILES_DIR}/common/ko_category.txt" --output_name "modules_category" --kegg_type "modules"

#python functional_enrichment/helpers/obtain_kegg_category.py --funit_transit_txt_file_path "${REF_FILES_DIR}/common/kegg_pathways.txt" --save_loc "${REF_FILES_DIR}/common" --ko_category_file "${REF_FILES_DIR}/common/ko_category.txt" --output_name "pathways_category" --kegg_type "pathways"


