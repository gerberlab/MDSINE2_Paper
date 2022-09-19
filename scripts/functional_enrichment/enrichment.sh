#!/bin/bash 
set -e 
source functional_enrichment/settings.sh 

#require_program picrust2 
echo "Running Picrust2. Saving output to ${PICRUST_OUTPUT_DIR}"

#echo "Step 1: Sequence Placement"
#place_seqs.py -s ${SEQUENCES_DIR} -o "${PICRUST_OUTPUT_DIR}/places_seqs.tre" -p 1 --intermediate working 

echo 

echo "Step 2: Hidden State Prediction"
echo "Hidden State Prediction for 16S copy"
#hsp.py -i 16S -t "${PICRUST_OUTPUT_DIR}/places_seqs.tre" -o "${PICRUST_OUTPUT_DIR}/marker_nsti_predicted.tsv.gz" -p 1 -n
#tar -zxvf "${PICRUST_OUTPUT_DIR}/marker_nsti_predicted.tsv.gz" 
echo 

echo "Hidden State Prediction for EC numbers"
#hsp.py -i EC -t "${PICRUST_OUTPUT_DIR}/places_seqs.tre" -o "${PICRUST_OUTPUT_DIR}/EC_predicted.tsv.gz" -p 1 
#tar -czvf "${PICRUST_OUTPUT_DIR}/EC_predicted.tsv.gz" "${PICRUST_OUTPUT_DIR}/EC_predicted.tsv" 
echo 

#echo "Hidden State Predictions for KO"
#hsp.py -i KO -t "${PICRUST_OUTPUT_DIR}/places_seqs.tre" -o "${PICRUST_OUTPUT_DIR}/KO_predicted.tsv.gz" -p 1 
#tar -czvf "${PICRUST_OUTPUT_DIR}/KO_predicted.tsv.gz" "${PICRUST_OUTPUT_DIR}/KO_predicted.tsv" 
echo 

echo "Step 3: Metagenome Prediction"
echo "Metagenome Prediction for EC"
metagenome_pipeline.py -i ${COUNTS_DIR}\
                       -m "${PICRUST_OUTPUT_DIR}/marker_nsti_predicted.tsv.gz" \
                       -f "${PICRUST_OUTPUT_DIR}/EC_predicted.tsv.gz" \
                       -o "${PICRUST_OUTPUT_DIR}/EC_metagenome_out"
echo

echo "Metagenome Prediction for KO"
metagenome_pipeline.py -i ${COUNTS_DIR}\
                       -m "${PICRUST_OUTPUT_DIR}/marker_nsti_predicted.tsv.gz" \
                       -f "${PICRUST_OUTPUT_DIR}/KO_predicted.tsv.gz" \
                       -o "${PICRUST_OUTPUT_DIR}/KO_metagenome_out"