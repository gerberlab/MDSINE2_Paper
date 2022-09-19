#!/bin/sh
set -e
source functional_enrichment/settings.sh

python functional_enrichment/helpers/obtain_ko_category.py --funit_ko_file "${REF_FILES_DIR}/common/kegg_modules.txt" --save_loc "${REF_FILES_DIR}/common" --save_name "modules_category.txt"
