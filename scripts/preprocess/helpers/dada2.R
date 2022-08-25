#if (!requireNamespace("BiocManager", quietly=TRUE))
#install.packages("BiocManager")
#BiocManager::install("dada2", force=TRUE)

library("optparse")
library(dada2); # we are using version 1.22
packageVersion("dada2")


option_list = list(
	make_option(c("-r", "--reads"), type="character", default=NULL, help="The directory containing the reads.", metavar="DIRNAME"),
    make_option(c("-o", "--outdir"), type="character", default=NULL, help="output file name [default= %default]", metavar="character"),
	make_option(c("-s", "--seed"), type="integer", default=100, help="The seed to use for random number generation.", metavar="INT"),
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# ==================================================== Initialization
set.seed(opt$seed)

path <- opt$reads
if (is.null(path)) {
	print_help(opt_parser)
	stop("The \"reads\" argument is required. See script usage (--help)")
}
list.files(path)

output_path <- opt$outdir
if (is.null(path)) {
	print_help(opt_parser)
	stop("The \"outdir\" argument is required. See script usage (--help)")
}
dir.create(output_path)

temp_path <- paste(output_path, '/tmp')
dir.create(temp_path)

rdp_train_path <- paste(output_path, '/taxa/rdp/rdp_train_set_18.fa.gz')
rdp_species_path <- paste(output_path, '/taxa/rdp/rdp_species_assignment_18.fa.gz')
silva_train_path <- paste(output_path, '/taxa/silva/silva_nr99_v138.fa.gz')
silva_species_path <- paste(output_path, '/taxa/silva/silva_species_assignment_v138.fa.gz')


# ==================================================== DATA PROCESSING

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

setDadaOpt(MAX_CONSIST=30)

errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE, nbases=1e8, pool = "pseudo")
errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE, nbases=1e8, pool = "pseudo")

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = "pseudo", selfConsist=TRUE) #19 consist steps
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = "pseudo", selfConsist=TRUE) #13 consist steps

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
