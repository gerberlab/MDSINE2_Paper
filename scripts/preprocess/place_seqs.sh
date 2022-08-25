#!/bin/bash
set -e
source preprocess/settings.sh

rdp_dir=$1
basename=$2
out_dir=$3
input_fasta=$4
start_idx=$5
end_idx=$6

require_variable 'rdp_dir' $rdp_dir
require_variable 'basename' $basename
require_variable 'out_dir' $out_dir
require_variable 'input_fasta' $input_fasta
require_variable 'start_idx' $start_idx
require_variable 'end_idx' $end_idx

require_program python
require_program hmmbuild
require_program hmmalign
require_program pplacer
require_program perl


# ===================== Main logic here
tmp_dir=${rdp_dir}/"__tmp"
mkdir -p $tmp_dir


echo "[*] Extracting target region ${start_idx}--${end_idx}."
subseq_file=${tmp_dir}/${basename}_${start_idx}_${end_idx}.fa
python preprocess/helpers/phylo_placement/trim_fasta.py \
--start $start_idx \
--end $end_idx \
-i ${rdp_dir}/${basename}_Aln.fa \
-o ${subseq_file}\


echo "[*] Creating HMM."
hmm_model=${tmp_dir}/${basename}_${start_idx}_${end_idx}.hmm
hmmbuild $hmm_model $subseq_file


echo "[*] Aligning sequences (${input_fasta})."
aln_output=${out_dir}/aligned_sequences.sto
hmmalign --trim --mapali ${subseq_file} -o $aln_output $hmm_model $input_fasta


echo "[*] Placing aligned sequences on the reference tree."
pplacer_output=${out_dir}/placement.jplace
pplacer --verbosity 1 -c . -o $pplacer_output $aln_output


echo "[*] Placing on XML tree."
xml_output=${tmp_dir}/xml_tree.xml
guppy tog -o $xml_output --xml $pplacer_output


echo "[*] Replace species id with taxonomy."
perl preprocess/helpers/phylo_placement/replace_id_taxaname.pl $xml_output ${rdp_dir}/${basename}_info.csv
xml_species_output=${tmp_dir}/xml_tree_speciesName.xml


echo "[*] Make the newick tree."
python preprocess/helpers/phylo_placement/phyloxml_to_newick.py -i $xml_output -o ${out_dir}/newick_tree_full_taxid.nhx
newick_output=${out_dir}/newick_tree_full_speciesName.nhx
python preprocess/helpers/phylo_placement/phyloxml_to_newick.py -i $xml_species_output -o $newick_output


echo "[*] Trim to only query reads."
python preprocess/helpers/phylo_placement/prune_tree.py \
-f $input_fasta \
-n $newick_output \
-o ${out_dir}/newick_tree_query_reads.nhx
