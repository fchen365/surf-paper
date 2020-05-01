## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## perform SURF analysis for one RBP
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## install R package
# install.packages("devtools")
devtools::install_github("fchen365/surf")
library(surf)

## ------- 1. parse -------
## parse ATR events from genome annotation
anno.file = "gencode.v24.annotation.filtered.gtf"
event <- parseEvent(anno.file, cores = 22) ## time-consuming, about 1 hour for filtered GTF
saveRDS(event, "encode_surf/GRCh38.v24.filtered.event.rds")

## ------- 2. drseq -------
sampleData.shrna <- data.frame(
  row.names = c('treated1', 'treated2', 'untreated1', 'untreated2'), ## sample names
  bam = paste0("rna-seq/bam/",1:4,".bam"), ## path to bam files
  condition = c('knockdown', 'knockdown', 'control', 'control'), 
  stringsAsFactors = F
)
drr <- drseq(event, sampleData.shrna, strandSpecific = 2) 

## ------- 3. faseq ------- 
sampleData.eclip = data.frame(
  row.names = c('IP1', 'IP2', 'SMI1'), 
  bam = paste0('clip-seq/bam/',5:7,'.bam'), 
  condition = c('IP', 'IP', 'SMI'),
  stringsAsFactors = F)
far <- faseq(drr, sampleData.eclip, strandSpecific = 2)

## ------- 4. daseq -------
## External data: TCGA vs GTEx
## rank transcripts (TPM)
exprMat <- readRDS('io/TcgaTargetGtex_rsem_isoform_tpm_laml_blood.rds')
rankings <- getRankings(exprMat)
## sample data
sampleData.external <- data.frame(
  condition = rep(c('TCGA', 'GTEx'), c(173, 337)),
  row.names = colnames(exprMat)) 
## differential activity (transcript)
dar <- daseq(getRankings(exprMat), sampleData.external, far)


saveRDS(dar, "encode_surf/RBP.NAME.results.rds")
## replace "RBP.NAME" with the actual RBP names
