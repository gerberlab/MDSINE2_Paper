#!/bin/bash
set -e
source runtime_benchmark/settings.sh

echo "[*] Compiling runtimes."

csv_path=${OUTPUT_DIR}/runtimes.csv
echo "MethodName,NumTaxa,Trial,Runtime" > $csv_path  # Header line


add_entry() {
	methodname=$1
	subdir=$2

	trial_dir=${OUTPUT_DIR}/taxa_top_${n_taxa}/trial_${trial}
	elapsed_time=$(cat ${trial_dir}/${subdir}/runtime.txt)
	echo "${methodname},${n_taxa},${trial},${elapsed_time}" >> $csv_path
}


for n_taxa in 10 25 50 100; do
	for (( trial = 1; trial < ${NUM_TRIALS}+1; trial++ )); do
		echo "Handling n_taxa: ${n_taxa}, trial: ${trial}"
		add_entry "MDSINE2(Negbin)" "mdsine2_negbin"
		add_entry "MDSINE2" "mdsine2"
		add_entry "cLV elastic-net" "clv/elastic-net"
		add_entry "gLV elastic-net" "glv/elastic-net"
		add_entry "gLV ridge" "glv/ridge"
		add_entry "gLV-ra elastic-net" "glv-ra/elastic-net"
		add_entry "gLV-ra ridge" "glv-ra/ridge"
		add_entry "lra" "lra/elastic-net"
done

echo "[*] Wrote results to ${csv_path}"
echo "[*] Done."
