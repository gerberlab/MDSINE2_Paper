#if (!requireNamespace("BiocManager", quietly=TRUE))
#install.packages("BiocManager")
#BiocManager::install("dada2", force=TRUE)

set.seed(100)
library(dada2); 
packageVersion("dada2")

path <- "C:/Users/sawal/Dropbox (Partners HealthCare)/all_input_data/healthy"


output_path <- '/output'
dir.create(output_path)


list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[10:11])
plotQualityProfile(fnRs[10:11])




