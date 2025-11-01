# Preprocessing of ASVs for MDSINE2

## Before you start!

This pipeline requires the following software: 
`biopython`, `ete3`, `six`, `guppy3`, `hmmer`, `pplacer`.

For example, install these using conda as shown below:

```bash
conda install -c conda-forge -c bioconda biopython ete3 six guppy3 hmmer pplacer
```

On Mac, it may be helpful to instead use a brew installation of pplacer:

```
brew install brewsci/bio/pplacer
```

## Step 1: Run DADA2

```bash
bash preprocess/run_dada2.sh
```


## Step 2: Run phylogenetic placement for ASVs (this also performs the necessary multiple alignment)

```bash
bash preprocess/phylo_placement_asv.sh healthy
```


## Step 3: Run agglomeration of ASVs into OTUs

```bash
bash preprocess/agglomerate_asvs.sh healthy
```


## Step 4: [Optional!] Run phylogenetic placement for OTUs

```bash
bash preprocess/phylo_placement_otus.sh healthy
```


## Step 5: Assign taxonomic labels to OTUs.

In general, RDP needs to be run separately on the input sequences, 
located in metadata_asv/prefiltered_asvs.fa.
In this repository, we've included pre-generated RDP outputs. Simply run the script below:

```bash
bash preprocess/assign_otu_taxonomy.sh healthy
```


## Step 6: Filter the OTUs.

```bash
bash preprocess/filter_otus.sh healthy
```


## Step 7: [Optional!] Visualize the OTUs' and their underlying ASVs' abundances.

```bash
bash preprocess/plot_otus_filtered.sh healthy
```

- TODO describe what is output by this pipeline.
