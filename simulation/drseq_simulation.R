## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## simulation study: DrSeq vs others
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressPackageStartupMessages({
  library(surf)
  library(cowplot)
  library(ggVennDiagram) 
})

meths <- c("DrSeq", "DEXSeq-L", "DEXSeq-C", "rMATS", "MAJIQ", "MAJIQ-L")
colrs <- c("#E69F00", '#66c2a5', "#56B4E9", "#984ea3", "#fccde5", "#fb8072") ## "#999999" (grey) "#0072B2" (dark blue) ## method colors


## ---- DrSeq ---- 
annotation.file <- 'Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf'

## parse events
event <- parseEvent(annotation.file, cores = 22) ## took 2.6 hours (unfiltered)

## drseq
sampleData = data.frame(
  bam = paste0("simulation/bam/Hs_sample",1:6,".bam"),
  condition = rep(c('I', 'II'), each = 3), 
  stringsAsFactors = F)
drr <- drseq(event, sampleData, cores = 6)

saveRDS(drr, "simulation/drseq/drseq.results.rds")

# ---- DEXSeq ----

countFiles = paste0("simulation/dexseq/Hs_sample",1:6,".txt")
sampleData.shrna = data.frame(
  row.names = 1:6,
  condition = rep(c('I', 'II'), each = 3), 
  stringsAsFactors = F)
flattenedFile = 'simulation/dexseq/GRCh37.71.dexseq.gff' 
dad <- DEXSeqDataSetFromHTSeq(
  countfiles = countFiles,
  sampleData = sampleData.shrna,
  ~ sample + exon + condition:exon,
  flattenedfile = flattenedFile
)
dxr <- DEXSeq(dad, BPPARAM = MulticoreParam(workers = 6))
saveRDS(dxr, "simulation/dexseq/dexseq.results.rds")


## ------ collect results ------

## ------ _ drseq (raw) ------ 
drr <- drseqResults(event)
res.drseq <- drr[c("gene_id", "transcript_id", "padj", "event_name")] %>% 
  as.data.frame %>% 
  mutate(test_id = rownames(.)) %>%
  filter(!is.na(padj))

table(event = drr$event_name[drr$padj < 0.05])
#   SE   RI A3SS A5SS  AFE  A5U  IAP  TAP
#  820  433  628  473 2575 1560 2377 2410

## ------ _ dexseq + adhoc ------ 
dxr = readRDS("simulation/dexseq/dexseq.results.rds")
res.dexseq <- dxr %>% as.data.frame %>% 
  transmute(test_id = rownames(.), gene_id = groupID, padj = padj) %>% 
  filter(!is.na(padj))

## adhoc
hit <- findOverlaps(drr$genomicData, dxr$genomicData)
padj <- split(dxr$padj[to(hit)], drr$event_id[from(hit)])
res.adhoc <- drr[c("gene_id", "transcript_id", "event_name")]
res.adhoc$test_id = rownames(res.adhoc)
res.adhoc$padj = padj[match(res.adhoc$test_id, names(padj))]
res.adhoc <- res.adhoc[!is.na(res.adhoc$padj),] %>% as.data.frame
## liberal
res.adhoc1 <- mutate(res.adhoc, padj = sapply(padj, min, na.rm = T))
## conservative
res.adhoc2 <- mutate(res.adhoc, padj = sapply(padj, quantile, .6, na.rm = T))


## ------ _ rmats ------ 
# list.files(pattern = "*JCEC.txt")
events <- c("SE", "RI", "A3SS", "A5SS", "MXE")
rmats <- foreach(e = events, .combine = rbind) %do% {
  file <- paste0("simulation/rmats/", e,".MATS.JCEC.txt")
  dat <- read.table(file, header = T, stringsAsFactors = F)
  cbind(dat[c("ID", "GeneID", "geneSymbol", "chr", "strand",
              "IncLevelDifference", "PValue", "FDR")], 
        event_name = e) 
}
res.rmats <- rmats %>% 
  transmute(test_id = paste(event_name, ID, sep = ":"), 
            gene_id = GeneID, 
            padj = FDR, 
            event_name = event_name)

## ------ _ majiq ------ 
## default
majiq <- read_tsv("simulation/majiq/I_II.deltapsi.2.tsv") ## default
nrow(majiq) ## 974
colSums(majiq[,c("A5SS", "A3SS", "ES")])
# A5SS A3SS   ES
#  177  198  850
res.majiq <- majiq %>% 
  transmute(
    test_id = `LSV ID`, 
    gene_id = `Gene ID`, 
    pvalue = `P(|dPSI|<=0.05) per LSV junction`, 
    padj = strsplit(pvalue, ";") %>% lapply(min) %>% 
      unlist() %>% as.numeric(),
    padj = 0,
    event_name = cbind(ifelse(ES, "SE", NA),
                       ifelse(A5SS, "A5SS", NA),
                       ifelse(A3SS, "A3SS", NA)) %>% 
      apply(1, function(x) 
        paste(na.omit(x), collapse = ","))
  )
dplyr::count(res.majiq, event_name)
#   event_name         n
# 1 ""                38
# 2 "A3SS"            52
# 3 "A5SS"            29
# 4 "A5SS,A3SS"        5
# 5 "SE"             596
# 6 "SE,A3SS"        111
# 7 "SE,A5SS"        113
# 8 "SE,A5SS,A3SS"    30

# 297/974 = 30.49%

sum(res.majiq$padj < .05) ## 791

## liberal version
majiq1 <- read_tsv("simulation/majiq/I_II.deltapsi.1.tsv") ## liberal
nrow(majiq) ## 2323
colSums(majiq1[,c("A5SS", "A3SS", "ES")])
# A5SS A3SS   ES
#  440  474 2044
res.majiq1 <- majiq1 %>% 
  transmute(
    test_id = `LSV ID`, 
    gene_id = `Gene ID`, 
    pvalue = `P(|dPSI|<=0.05) per LSV junction`, 
    padj = strsplit(pvalue, ";") %>% lapply(min) %>% 
      unlist() %>% as.numeric(),
    padj = 0,
    event_name = cbind(ifelse(ES, "SE", NA),
                       ifelse(A5SS, "A5SS", NA),
                       ifelse(A3SS, "A3SS", NA)) %>% 
      apply(1, function(x) 
        paste(na.omit(x), collapse = ","))
  )
dplyr::count(res.majiq1, event_name)
#   event_name         n
# 1 ""                86
# 2 "A3SS"            97
# 3 "A5SS"            78
# 4 "A5SS,A3SS"       18
# 5 "SE"            1414
# 6 "SE,A3SS"        286
# 7 "SE,A5SS"        271
# 8 "SE,A5SS,A3SS"    73

# 734 / 2323 = 31.60%

sum(res.majiq1$padj < .05) ## 2039

## ------ _ all together ------ 
res <- list(DrSeq = res.drseq, 
            "DEXSeq-L" = res.adhoc1, 
            "DEXSeq-C" = res.adhoc2, 
            rMATS = res.rmats, 
            MAJIQ = res.majiq, 
            "MAJIQ-L" = res.majiq1) %>% 
  lapply(magrittr::extract, c("test_id", "gene_id", "event_name", "padj")) %>% 
  bind_rows(.id = "method") %>% 
  mutate(method = factor(method, meths)) 

## ------ stats: event type ------ 

## # of significant events
res %>% 
  filter(!grepl("MAJIQ", method)) %>%
  mutate(event_name = factor(event_name, c(surf.events, "MXE"))) %>%
  group_by(method, event_name) %>% 
  summarise(n_detected = sum(padj < .05, na.rm = T)) %>%
  ggplot(aes(event_name, n_detected, fill = method)) +
  geom_bar(stat = "identity",color = "grey60", size = .2, alpha = .9, 
           position = position_dodge(preserve = "single")) +
  scale_y_continuous(trans = "sqrt", breaks = c(.1,1,2.5,5)*1e3) + 
  labs(y = "# of detected events") +
  scale_fill_manual(values = setNames(colrs, meths)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) + 
  ggsave("bar_drseq_simu_detected_event.pdf", width = 7, height = 3)

#            SE  RI A3SS A5SS  AFE  A5U  IAP  TAP MXE
# DrSeq     820 433  628  473 2575 1560 2377 2410   0
# DEXSeq-L  937 495  642  530 4652 1745 5712 2730   0
# DEXSeq-C  538 250  431  384  782  967  827 1388   0
# rMATS    1662 462  354  250    0    0    0    0  74


## venn: MAJIQ detected event
g.venn <- filter(res, grepl("MAJIQ", method)) %>% 
  transmute(method = method, 
            SE = grepl("SE", event_name), 
            A3SS = grepl("A3SS", event_name),
            A5SS = grepl("A5SS", event_name)) %>% 
  group_split(method, keep = F) %>% 
  lapply(function(x) {
    lapply(x, which) %>% 
      ggVennDiagram(label = "count", show.legend = F) + 
      scale_fill_distiller(trans = "log", palette = "GnBu", direction = 1)
  })
plot_grid(g.venn[[1]], g.venn[[2]], ncol = 2, labels = "auto") + 
  ggsave("venn_majiq_detected_event.pdf", width = 7, height = 3)


## ------ TPR & FDR ------
## ------ _ truth ------ 
tpm <- sapply(seq_len(6), function(i) {
  file = paste0("simulation/truth/Hs_sample",i,"_TPM.txt")
  tpm <- read.table(file, header = T, stringsAsFactors = F)
  setNames(tpm$TPM, tpm$transcript_id)
})
truth <- read.table("simulation/truth/Hs_truth.txt", 
                    header = T, stringsAsFactors = F) %>%
  mutate(tpm = rowMeans(tpm)) %>% 
  group_by(gene_id) %>%
  mutate(rank = rank(-tpm),
         tx_ds_status = gene_ds_status * (rank <= 2), 
         IsoPct = tpm / sum(tpm)) 


## this table is used by the following two functions 
truth_gene <- truth %>% 
  filter(!is.nan(IsoPct)) %>% ## filter out genes with sum(tpm)==0
  # filter(tx_ds_status == 1) %>% ## filter by selected tx
  group_by(gene_id, gene_ds_status) %>% 
  summarise(diffIsoPct = max(IsoPct[rank <= 2]) - min(IsoPct[rank <= 2])) %>% 
  ungroup() %>% 
  mutate(diffIsoPct = diffIsoPct * gene_ds_status, 
         diffIsoPct = ifelse(diffIsoPct, diffIsoPct, runif(length(diffIsoPct), 0, max(diffIsoPct))))

## FDR calculator (individual tests as units)
fdr <- function(gene_id, padj, alpha = 0.05, 
                minDiffIsoPct = 0, maxDiffIsoPct = 1){
  padj[is.na(padj)] = 1 ## treat NA as 1
  sub_truth_gene <- truth_gene %>% 
    filter(diffIsoPct >= minDiffIsoPct, 
           diffIsoPct <= maxDiffIsoPct)
  tested_gene_id <- sub_truth_gene$gene_id
  truth_gene_id <- sub_truth_gene$gene_id[sub_truth_gene$gene_ds_status == 1]
  p <- (gene_id %in% tested_gene_id) & (padj < alpha)
  tp <- (gene_id %in% truth_gene_id) & (padj < alpha)
  return((sum(p) - sum(tp)) / sum(p))
}

## TPR calculator (individual gene as units)
tpr <- function(gene_id, padj, alpha = 0.05, 
                minDiffIsoPct = 0, maxDiffIsoPct = 1){
  padj[is.na(padj)] = 1 ## change NA as 1
  sub_truth_gene <- truth_gene %>% 
    filter(diffIsoPct >= minDiffIsoPct, 
           diffIsoPct <= maxDiffIsoPct)
  tested_gene_id <- sub_truth_gene$gene_id
  truth_gene_id <- sub_truth_gene$gene_id[sub_truth_gene$gene_ds_status == 1]
  tp <- (gene_id %in% truth_gene_id) & (padj < alpha)
  n_tp <- length(unique(gene_id[tp]))
  n_t <- length(unique(truth_gene_id))
  return(n_tp / n_t)
}


## ------ _ overall ------ 
g1 <- lapply(c(0.1, 0.05, 0.01), function(alpha) {
  res %>% 
    filter(!grepl("MAJIQ", method)) %>%
    group_by(method) %>% 
    summarise(FDR = fdr(gene_id, padj, alpha = alpha), 
              TPR = tpr(gene_id, padj, alpha = alpha), 
              alpha = alpha)
}) %>% bind_rows() %>% 
  ggplot(aes(FDR, TPR, color = method)) + 
  geom_vline(xintercept = c(0.01, 0.05, 0.1), 
             color = "grey50", linetype = 2, alpha = .8) +
  geom_point(aes(shape = method), alpha = .9) + 
  geom_line(alpha = 0.8) + 
  scale_color_manual(values = setNames(colrs, meths)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  ggsave("scatter_drseq_simu_TPRvsFDR_overall.pdf", width = 4.2, height = 3)

## ------ _ by event type ------ 
## stratefied 
g2 <- lapply(c(0.1, 0.05, 0.01), function(alpha) {
  res %>% 
    filter(!grepl("MAJIQ", method)) %>%
    mutate(event = factor(event_name, c(surf.events, "MXE"))) %>%
    group_by(method, event) %>% 
    summarise(FDR = fdr(gene_id, padj, alpha = alpha), 
              TPR = tpr(gene_id, padj, alpha = alpha), 
              alpha = alpha)
}) %>% bind_rows() %>% 
  ggplot(aes(FDR, TPR, color = method)) + 
  geom_vline(xintercept = c(0.01, 0.05, 0.1), 
             color = "grey50", linetype = 2, alpha = .8) + 
  geom_point(aes(shape = method), alpha = .9) + 
  geom_line(alpha = 0.8) + 
  scale_color_manual(values = setNames(colrs, meths)) +
  facet_wrap(~ event, ncol = 3) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  ggsave("scatter_drseq_simu_TPRvsFDR_stratefied_event.pdf", width = 7, height = 6)


## ------ _ by diff IsoPct ------ 
g3 <- lapply(c(0.1, 0.05, 0.01), function(alpha) {
  lapply(c("small", "medium", "large"), function(delta) {
    minDiffIsoPct <- switch(delta,
                            small = 0,
                            medium = .05,
                            large = .1)
    maxDiffIsoPct <- switch(delta,
                            small = .05,
                            medium = .1,
                            large = 1)
    res %>% 
      group_by(method) %>% 
      summarise(FDR = fdr(gene_id, padj, alpha = alpha, 
                          minDiffIsoPct = minDiffIsoPct, 
                          maxDiffIsoPct = maxDiffIsoPct), 
                TPR = tpr(gene_id, padj, alpha = alpha, 
                          minDiffIsoPct = minDiffIsoPct, 
                          maxDiffIsoPct = maxDiffIsoPct), 
                diffIsoPct = paste0("(", minDiffIsoPct, ", ", 
                                    maxDiffIsoPct, "]"))
  }) %>% bind_rows() 
}) %>% bind_rows() %>% 
  ggplot(aes(FDR, TPR, color = method)) + 
  geom_vline(xintercept = c(0.01, 0.05, 0.1), 
             color = "grey50", linetype = 2, alpha = .8) + 
  geom_point(aes(shape = method), alpha = .9) + 
  geom_line(alpha = 0.8) + 
  scale_x_continuous(breaks = c(0,.3,.6,.9)) + 
  scale_color_manual(values = setNames(colrs, meths)) +
  facet_wrap(~ diffIsoPct, ncol = 3, 
             labeller = label_bquote(Delta*"IsoPct"%in%.(diffIsoPct))) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  ggsave("scatter_drseq_simu_TPRvsFDR_stratefied_IsoPct.pdf", 
         width = 8, height = 2.4)

