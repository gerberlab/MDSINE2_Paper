

```bash
conda create -n phyloplacement_py3
conda activate phyloplacement_py3
pip install biopython
pip install ete3
pip install six
pip install guppy3
conda install -c bioconda hmmer
conda install -c bioconda pplacer #on Mac: brew install brewsci/bio/pplacer
```




Place ASVs on pre-aligned RDP high quality 16S reads

```bash
python scripts/place_seqs.py \
    --v4-region-start 1045 \
    --v4-region-end 1374 \
    --refpkg RDP-11-5_TS_Processed.refpkg \
    --query-reads ../sequence_analysis/output/sequences.fa \
    --output-folder output_ASVs \
    --temp-folder tmp_ASVs/
```

Create alignment file for ASVs

```bash
python scripts/sto_to_fasta.py -i output_ASVs/placed_sequences_on_v4_region.sto -o aligned_asvs/aligned_asvs.fa
```

Agglomerate ASVs
