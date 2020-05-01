## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## prepare TCGA and GTEx transcriptome data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ------ import from xena (TCGA+TARGET+GTEx) ------
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

## filter LAML + whole blood 
pheno = read.delim('xena/TcgaTargetGtex_phenotype.txt', 
                   stringsAsFactors = F)
pheno_laml = pheno %>% filter(
  (X_study == 'GTEX' & 
     primary.disease.or.tissue == "Whole Blood") | 
    (X_study == 'TCGA' & 
       primary.disease.or.tissue == "Acute Myeloid Leukemia"))


## isoform TPM
tbl = fread('xena/TcgaTargetGtex_rsem_isoform_tpm', data.table = F)
rownames(tbl) = as.character(tbl$sample)
tbl = tbl[,-1]
tbl_laml = tbl[,pheno_laml$sample]
saveRDS(tbl_laml, 'io/TcgaTargetGtex_rsem_isoform_tpm_laml_blood.rds')


## gene TPM
tbl = fread('xena/TcgaTargetGtex_rsem_gene_tpm', data.table = F)
rownames(tbl) = as.character(tbl$sample)
tbl = tbl[,-1]
tbl_laml = tbl[,pheno_laml$sample]
saveRDS(tbl_laml, 'io/TcgaTargetGtex_rsem_gene_tpm_laml_blood.rds')


## isoform percentage
tbl = fread('xena/TcgaTargetGtex_rsem_isopct', data.table = F)
rownames(tbl) = as.character(tbl$sample)
tbl = tbl[,-1]
tbl_laml = tbl[,pheno_laml$sample]
saveRDS(tbl_laml, 'io/TcgaTargetGtex_rsem_isopct_laml_blood.rds')
