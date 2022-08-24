# Preprocessing of ASVs for MDSINE2

## Step 1: Run DADA2

- Specify how to obtain the reads. Should we host it somewhere on the web?

- TODO describe DADA2 script/how to use.

## Step 2: Run phylogenetic placement for ASVs (this also performs the necessary multiple alignment)

- todo create a script which does this for the user.
```bash
python scripts/place_seqs.py \
    --v4-region-start 1045 \
    --v4-region-end 1374 \
    --refpkg RDP-11-5_TS_Processed.refpkg \
    --query-reads ../sequence_analysis/output/sequences.fa \
    --output-folder output_ASVs \
    --temp-folder tmp_ASVs/
```

```bash
python scripts/sto_to_fasta.py -i output_ASVs/placed_sequences_on_v4_region.sto -o aligned_asvs/aligned_asvs.fa
```

## Step 3: Run agglomeration of ASVs into OTUs

```bash
bash preprocess/agglomerate_asvs.sh
```

## Step 4: [Optional!] Run phylogenetic placement for OTUs

- todo create a script which does this for the user.

```bash
python scripts/place_seqs.py \
    --v4-region-start 1045 \
    --v4-region-end 1374 \
    --refpkg RDP-11-5_TS_Processed.refpkg \
    --query-reads ../datasets/gibson/preprocessed/gibson_healthy_agg.fa \
    --output-folder output_OTUs \
    --temp-folder tmp_OTUs/
```

## Step 5: Assign taxonomic labels to OTUs.

- todo specify that RDP needs to be run and outputs downloaded to the correct folder.
- does RDP have a web API? (REST/AJAX?)

```bash
bash preprocess/assign_otu_taxonomy.sh
```

## Step 6: Filter the OTUs.

```bash
bash preprocess/filter_otus.sh
```

# Step 7: [Optional!] Visualize the OTUs' and their underlying ASVs' abundances.

```bash
bash preprocess/plot_otus_filtered.sh
```

- TODO describe what is output by this pipeline.
