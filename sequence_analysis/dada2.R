#if (!requireNamespace("BiocManager", quietly=TRUE))
#install.packages("BiocManager")
#BiocManager::install("dada2", force=TRUE)

set.seed(100)
library(dada2); # we are using version 1.22
packageVersion("dada2")

path <- './raw_reads/healthy'
output_path <- './output'
temp_path <- './tmp'

rdp_train_path <- './taxa/rdp/rdp_train_set_18.fa.gz'
rdp_species_path <- './taxa/rdp/rdp_species_assignment_18.fa.gz'

silva_train_path <- './taxa/silva/silva_nr99_v138.1_train_set.fa.gz'
silva_species_path <- './taxa/silva/silva_species_assignment_v138.1.fa.gz'


dir.create(temp_path)
dir.create(output_path)

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[11:12])
plotQualityProfile(fnRs[11:12])

filtFs <- file.path(temp_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(temp_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(245,150), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

head(out)

setDadaOpt(PSEUDO_PREVALENCE=4)

errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE, nbases=5e8, pool = "pseudo")
errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE, nbases=5e8, pool = "pseudo")

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = "pseudo")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

# Remove sequeunces that are much longer or shorter than expected due to non-specific priming.
# seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)

# Create a table that tracks reads throughout the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, file.path(output_path, "track.tsv"), sep='\t', col.names=NA)

# Writing files
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# Write a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta,  file.path(output_path, "sequences.fa"))


asv_headers <- sub(">", "", asv_headers)
asv_tab <- t(seqtab.nochim)
rownames(asv_tab) <- asv_headers
write.table(asv_tab, file.path(output_path, "counts.tsv"), sep="\t", col.names=NA)

taxa_rdp <- assignTaxonomy(seqtab.nochim, rdp_train_path, multithread=TRUE)
taxa_rdp <- addSpecies(taxa_rdp, rdp_species_path, allowMultiple=2)
taxa_rdp <- cbind(rownames(taxa_rdp), taxa_rdp)

rownames(taxa_rdp) <- asv_headers
colnames(taxa_rdp) <- c("sequence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
write.table(taxa_rdp, file.path(output_path, "rdp_species.tsv"), sep="\t", col.names=NA)


taxa_silva <- assignTaxonomy(seqtab.nochim, silva_train_path, multithread=TRUE)
taxa_silva <- addSpecies(taxa_silva, silva_species_path, allowMultiple=2)
taxa_silva <- cbind(rownames(taxa_silva), taxa_silva)
write.table(taxa_silva, file.path(output_path, "silva_species.tsv"), sep="\t", col.names=NA)
