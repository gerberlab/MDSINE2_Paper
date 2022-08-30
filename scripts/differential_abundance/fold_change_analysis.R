rm(list=ls())

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("apeglm")
library("DESeq2") #v1.32.0 was used for this paper


taxa_rank<-c('phylum', 'of')# 'otu', 'of','phylum'
cohort_name<-c('healthy', 'uc')
folder<-'out/'

for (cohort in cohort_name)
{
  for (taxa in taxa_rank)
  {
  name<-paste(cohort,taxa, sep='_')
  print(name)
  print(paste(folder,name,'_counts.csv', sep=''))
  
  cts <- read.csv(paste(folder,name,'_counts.csv', sep=''), row.names = 1)
  coldata <- read.csv(paste(folder,name,'_meta.csv', sep=''))

  coldata$window<- factor(coldata$window)

  dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~window)
                          
  akeep <- rowSums(counts(dds)) >= 100 
  dds <- dds[akeep,]
  dds<- DESeq(dds)
  
  contrast <- c("window", "1", "0")
  res<- results(dds, contrast=contrast)
  resOrdered <- res[order(res$pvalue),]
  write.csv(as.data.frame(resOrdered), file=paste(folder,name,'_window_1.csv', sep=''))

  contrast <- c("window", "2", "1.5")
  res<- results(dds, contrast=contrast)
  resOrdered <- res[order(res$pvalue),]
  write.csv(as.data.frame(resOrdered), file=paste(folder,name,'_window_2.csv', sep=''))


  contrast <- c("window", "3", "2.5")
  res<- results(dds, contrast=contrast)
  resOrdered <- res[order(res$pvalue),]
  write.csv(as.data.frame(resOrdered), file=paste(folder,name,'_window_3.csv', sep=''))

  }
}

?results

for (taxa in taxa_rank)
{
  name<-paste('ss',taxa, sep='_')
  print(name)
  print(paste(folder,name,'_counts.csv', sep=''))
  
  cts <- read.csv(paste(folder,name,'_counts.csv', sep=''), row.names = 1)
  coldata <- read.csv(paste(folder,name,'_meta.csv', sep=''))
  
  
  coldata$cohort<- factor(coldata$cohort)
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~cohort)
  
  akeep <- rowSums(counts(dds)) >= 100 
  dds <- dds[akeep,]
  dds<- DESeq(dds)
  
  contrast <- c("cohort", "uc", "healthy")
  res<- results(dds, contrast=contrast)
  resOrdered <- res[order(res$pvalue),]
  write.csv(as.data.frame(resOrdered), file=paste(folder,name,'.csv', sep=''))
}




