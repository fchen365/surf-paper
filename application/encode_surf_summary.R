## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## summarise SURF analysis on ENCODE data 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressPackageStartupMessages({
  library(surf)
  library(scales)
  library(ggtree)
  library(ggpubr)
  library(ggrepel)
  library(cowplot)
  library(pheatmap)
  library(reshape2)
  library(universalmotif)
})

## 104 RBP names 
targets <- read.delim("data/RBP.txt", stringsAsFactors = F)$RBP 
functions <- read.delim("data/RBP.txt", row.names = "RBP")
functions$Function = factor(functions$Function, unique(functions$Function))
selected.events <- c("SE", "RI", "A3SS", "A5SS")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ ATR events (module 0) ------- 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
anno_event <- readRDS("encode_surf/GRCh38.v24.filtered.event.rds")

## ------ _ event proportion ------
nrow(anno_event) ## 199877

table(type = anno_event$event_name)
#    SE    RI  A3SS  A5SS   AFE   A5U   IAP   TAP
# 21599 12208 19076 16447 32138 34971 29748 33690
(32138 + 34971 + 29748 + 33690) / 199877
21599 / 199877 ## SE, 0.1080615
34971 / 199877 ## A5U, 0.1749626
33690 / 199877 ## TAP, 0.1685537

## bar
table(type = anno_event$event_name) %>% 
  data.frame() %>% 
  mutate(type = factor(type, surf.events)) %>%
  ggplot(aes(x = "%", Freq, fill = type)) + 
  geom_bar(stat = "identity", position = "fill", show.legend = F) +
  labs(x = "", y = "") +
  scale_fill_manual(values = setNames(surf.colors, surf.events)) + 
  scale_y_continuous(labels = percent) +
  coord_flip() + 
  theme_pubr() + 
  theme(axis.line = element_blank(), 
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggsave("bar_event_dist.pdf", width = 5, height = .7)


## ------ _ hist: #events per gene ------
table(gene = anno_event$gene_id, 
      event = anno_event$event_name) %>% 
  reshape2::melt() %>% 
  mutate(count = cut(value, breaks = c(-Inf, 0,1,2,3,4,5,10,20,Inf), 
                     labels = c(0:5, "[6-10]", "[11-20]", ">20"))) %>%
  ggplot(aes(count, fill = count)) + 
  geom_bar(stat = "count", color = "grey50", size = 0.2, alpha = .9) +
  facet_wrap(~ event, ncol = 4) + 
  labs(x = "# of ATR events", y = "# of genes") +
  scale_fill_brewer(palette = "GnBu") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  ggsave("bar_n_event_per_gene.pdf", width = 8, height = 4)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ DrSeq (module 1) ------- 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------- _ mean-dispersion functions ------- 
g1 <- readRDS("rna-seq/drseq/SRSF1.results.rds") %>% plotDispFunc() + 
  scale_y_continuous(limits = c(1.5e-4, 5), breaks = c(0.001, 0.01, 0.1, 1), trans = "log10")
g2 <- readRDS("rna-seq/drseq/AQR.results.rds") %>% plotDispFunc() + 
  scale_y_continuous(limits = c(1.5e-4, 5), breaks = c(0.001, 0.01, 0.1, 1), trans = "log10") 
g3 <- readRDS("rna-seq/drseq/HNRNPC.results.rds") %>% plotDispFunc() + 
  scale_y_continuous(limits = c(1.5e-4, 5), breaks = c(0.001, 0.01, 0.1, 1), trans = "log10") 
g4 <- readRDS("rna-seq/drseq/SF3B4.results.rds") %>% plotDispFunc() + 
  scale_y_continuous(limits = c(1.5e-4, 5), breaks = c(0.001, 0.01, 0.1, 1), trans = "log10") 
g <- ggarrange(g1, g2, g3, g4, labels = "auto", align = "hv", nrow = 2, ncol = 2) 
ggsave("line_mean_dispersion.pdf", width = 7, height = 7)


## ------- _ volcano plot ------- 
g = readRDS("rna-seq/drseq/AQR.results.rds") %>% 
  volcano.plot(remove.portion.grey =.9) +
  ggsave("volcano_AQR.pdf", width = 8, height = 4)
g = readRDS("rna-seq/drseq/SRSF1.results.rds") %>% 
  volcano.plot(remove.portion.grey =.9) + 
  ggsave("volcano_SRSF1.pdf", g, width = 8, height = 4)
g = readRDS("rna-seq/drseq/CPSF6.results.rds") %>% 
  volcano.plot(remove.portion.grey =.9) + 
  ggsave("volcano_CPSF6.pdf", g, width = 8, height = 4)
g = readRDS("rna-seq/drseq/SF3B4.results.rds") %>% 
  volcano.plot(y.limits = c(0, 20), remove.portion.grey =.9) + 
  ggsave("volcano_SF3B4.pdf", g, width = 8, height = 4)


## ------- _ dist of raw p-values ------- 
g <- c("SRSF1", "AQR", "HNRNPC", "SF3B4") %>% 
  lapply(function(target) {
    drr <- readRDS(paste0("rna-seq/drseq/", target, ".results.rds"))
    data.frame(p.value = drr$pvalue) %>% 
      ggplot(aes(p.value, fill = cut(p.value, 9))) +
      geom_histogram(bins = 50, alpha = .9, show.legend = FALSE) +
      labs(x = paste0("p-value"), y = "# of observations (thousand)") + 
      scale_fill_brewer(palette = "GnBu", direction = -1) + 
      scale_y_continuous(labels = function(x) x*1e-3) + 
      theme_bw()
  })
ggarrange(plotlist = g, ncol = 2, nrow = 2, 
          labels = "auto", common.legend = T) + 
  ggsave("hist_drseq_pval.pdf", width = 7, height = 7)


## ------- _ collect DR event_id ------- 
registerDoParallel(16)
event_ids_raw = foreach(target = targets) %dopar% {
  sdat = readRDS(paste0("encode_surf/", target, ".data.rds"))
  na.omit(sdat$event_id[sdat$padj < .05])
}
event_ids_filtered = foreach(target = targets) %dopar% {
  sdat = readRDS(paste0("encode_surf/", target, ".data.rds"))
  na.omit(sdat$event_id[sdat$padj < .05 & sdat$adjMean > 0.05])
}
event_ids_surf = foreach(target = targets) %dopar% {
  sdat = readRDS(paste0("encode_surf/", target, ".data.rds"))
  sdat$event_id[sdat$group != "no change" & sdat$included]
}
stopImplicitCluster()
names(event_ids_raw) = 
  names(event_ids_filtered) = 
  names(event_ids_surf) = targets
length(unique(unlist(event_ids_raw))) ## total #{DR events}=159673
length(unique(unlist(event_ids_filtered))) ## 130145
length(unique(unlist(event_ids_surf))) ## 80402

save(event_ids_raw, event_ids_filtered, event_ids_surf, 
     file = "encode_surf/drseq.data.all.RData")

## ------- __ event dist (overall) ------- 
load("encode_surf/drseq.data.all.RData")
# event_ids <- event_ids_raw
# event_ids <- event_ids_filtered
event_ids <- event_ids_surf
## aggregate 104 RBP
anno_event <- readRDS("encode_surf/GRCh38.v24.filtered.event.rds")
event_dr = anno_event[unique(unlist(event_ids, use.names = F)),]
table(type = event_dr$event_name) %>% 
  data.frame %>% 
  mutate(type = factor(type, surf.events)) %>%
  ggplot(aes("%", Freq, fill = type)) + 
  geom_bar(stat = "identity", position = "fill", show.legend = F) +
  labs(x = "", y = "") +
  scale_fill_manual(values = setNames(surf.colors, surf.events)) + 
  scale_y_continuous(labels = percent) +
  coord_flip() + 
  theme_pubr() + 
  theme(axis.line = element_blank(), 
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggsave("bar_event_dr_dist.pdf", width = 5, height = .7)


## ------- __ event dist (each RBP) ------- 

## how many increase/decrease on average
registerDoParallel(16)
drCountByGroup = foreach(target = targets) %dopar% {
  paste0("encode_surf/", target, ".data.rds") %>% 
    readRDS() %>% 
    .$group[.$group != "no change" & .$included] %>% 
    table()
} %>% 
  setNames(targets) %>%
  bind_rows(.id = "factor") %>% 
stopImplicitCluster()
colMeans(drCountByGroup)

## dr event count for each RBP
registerDoParallel(16)
drCount = foreach(target = targets) %dopar% {
  paste0("encode_surf/", target, ".data.rds") %>% 
    readRDS() %>% 
    .[.$group != "no change" & .$included, ## this is "event_ids_surf"
      c("event_id", "event_name", "logFoldChange")] %>% 
    data.frame() %>% 
    dplyr::rename("event" = "event_name") %>%
    mutate(LFC = cut(abs(logFoldChange), c(0, .1, .5, 1, Inf)), 
           LFC = recode(LFC, "(1,Inf]" = ">1")) %>% 
    group_by(event, LFC) %>% 
    summarise(count = n()) %>%
    ungroup() %>% 
    mutate(percent = sum(count[event %in% selected.events]) / sum(count))
} %>% 
  setNames(targets) %>%
  bind_rows(.id = "factor") %>% 
  arrange(percent, factor) %>% 
  filter(percent > 0, percent < 1) %>%
  mutate(factor = factor(factor, unique(factor)), 
         group = functions[as.character(factor), "Function"])
stopImplicitCluster()

## bar: event type proportion
g1 <- drCount %>% 
  group_by(factor, event, group) %>%
  summarise(count = sum(count)) %>%
  ggplot(aes(factor, count, fill = event)) +
  geom_bar(stat = "identity", position = "fill", 
           color = "grey90", size = 1e-3,
           width = .7, show.legend = F) + 
  geom_hline(yintercept = c(0,1), color = "grey90") +
  labs(x = "", y = "") + 
  scale_fill_manual(values = setNames(surf.colors, surf.events)) + 
  scale_y_continuous(labels = percent) +
  facet_grid(rows = vars(group), switch = "y", scales = "free_y", space = "free_y") + 
  coord_flip() + 
  theme_pubr() + 
  theme(panel.spacing.y = unit(.1, "lines"),
        strip.placement = "outside", 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 5), 
        axis.text = element_blank(), 
        axis.ticks.x = element_blank()) + 
  ggsave("bar_perc_dr_event.pdf", width = 5, height = 9)

## bar: number/count of DR events
g2 <- drCount %>% 
  group_by(factor, LFC, group) %>%
  summarise(count = sum(count)) %>%
  group_by(factor) %>% 
  arrange(LFC) %>%
  mutate(count.sqrt = sqrt(cumsum(count)), 
         count.sqrt = c(count.sqrt[1], diff(count.sqrt))) %>%
  ggplot(aes(factor, count.sqrt, fill = LFC)) +
  geom_bar(stat = "identity", color = "grey70", size = .1,
           width = .7, show.legend = T, 
           position = position_stack(reverse = T)) + 
  labs(y = "# of differential ATR events") + 
  scale_fill_brewer(palette = "GnBu") + 
  scale_y_continuous(labels = function(x) x^2) + 
  facet_grid(rows = vars(group), scales = "free_y", space = "free_y") + 
  coord_flip() + 
  theme_pubr() + 
  theme(panel.spacing.y = unit(.1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank()) + 
  ggsave("bar_dr_count_lfc.pdf", width = 3, height = 9)
g0 <- get_legend(
  g2 + theme_bw() + guides(
    title.position = "top",
    # direction = "horizontal",
    frame.colour = "grey70"
  ) + labs(fill = bquote("log"[2]~"fold change")) 
)
ggsave("guide_dr_lfc.pdf", g0, width = 2, height = 2)
g2 <- g2 + theme(legend.position = "none")

ggarrange(g1, g2, widths = c(3.5, 1.5), align = "h") + 
  ggsave("combo_event_dr_rbp.pdf", width = 7, height = 9)

## ------- _ compare w. rMATS ------- 
## compare # of detected event
## rMATS
rmats <- lapply(targets, function(target) {
  folder <- list.files("rna-seq/rmats", paste0(target,"-.*-K562"))
  if (!length(folder)) {
    message("No data for ", target, ", skipped.")
    return(data.frame())
  } 
  lapply(selected.events, function(e) {
    paste0("rna-seq/rmats/", folder, "/MATS_output/",
           e,".MATS.JunctionCountOnly.txt") %>%
      read.table(header = T, stringsAsFactors = F) %>% 
      filter(FDR < 0.05) %>% 
      transmute(test_id = paste0(target, "-", e, ":", ID), 
                gene_id = GeneID, 
                factor = target,  
                event = factor(e, selected.events))
  }) %>% bind_rows()
}) %>% bind_rows() %>% 
  group_by(factor, event) %>% 
  summarise(n_event = length(unique(test_id)), 
            n_gene = length(unique(gene_id)))

## DrSeq
anno_event <- readRDS("encode_surf/GRCh38.v24.filtered.event.rds")
load("encode_surf/drseq.data.all.RData")
drseq <- lapply(event_ids_surf, function(id) {
  anno_event[id, c("event_name", "event_id", "gene_id")] %>% 
    setNames(c("event", "test_id", "gene_id")) %>% 
    data.frame() %>% 
    filter(event %in% selected.events) %>% 
    mutate(event = factor(event, selected.events))
}) %>% 
  bind_rows(.id = "factor") %>% 
  group_by(factor, event) %>% 
  summarise(n_event = length(unique(test_id)), 
            n_gene = length(unique(gene_id)))

## scatter: # of detected events
dat <- bind_rows("DrSeq" = drseq, 
          "rMATS" = rmats, 
          .id = "method") %>% 
  reshape2::dcast(factor + event ~ method, value.var = "n_event")

## scatter: # of detected event
replace(dat, is.na(dat), .8) %>% ## this is for plotting 0 value
  ggplot(aes(rMATS, DrSeq, color = event)) + 
  geom_point(alpha = .6, show.legend = F) + 
  geom_rug(size = .2, alpha = .7, sides = "rt", show.legend = F) + 
  geom_abline(slope = 1, color = "grey60", linetype = 2) + 
  scale_x_continuous("# of events detected by rMATS", 
                     labels = scales::comma, trans = "log") +
  scale_y_continuous("# of events detected by DrSeq", 
                     labels = scales::comma, trans = "log") +
  scale_color_manual(values = setNames(surf.colors, surf.events)) + 
  facet_wrap(~ event, nrow = 2, scales = "free") + 
  theme_bw() + 
  ggsave("scatter_drseq_rmats_event.pdf", width = 6, height = 5.4)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ FASeq (module 2) ------- 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha = .05 ## FDR cut-off

# collect faseq results
registerDoParallel(21)
fit = foreach(target = targets) %dopar% {
  paste0('io/encode_surf/', target, '.results.rds') %>%
    readRDS() %>%
    faseqResults() %>%
    data.frame()
} %>%
  setNames(targets) %>%
  bind_rows(.id = "factor") %>%
  mutate(factor = as.factor(factor))
stopImplicitCluster()

## ------ _ FA plot ------
readRDS("encode_surf/AQR.results.rds") %>% 
  fa.plot(plot.event = c("SE", "AFE", "A5U")) + 
  ggsave("fa.plot_AQR.pdf", width = 8, height = 5)
readRDS("encode_surf/AQR.results.rds") %>%
  fa.plot(plot.event = c("SE")) +
  ggsave("fa.plot_AQR_SE.pdf", width = 5, height = 5)
readRDS("encode_surf/CPSF6.results.rds") %>%
  fa.plot(plot.event = "TAP") +
  ggsave("fa.plot_CPSF6.pdf", width = 4, height = 5)
readRDS("encode_surf/GEMIN5.results.rds") %>% 
  fa.plot(plot.event = c("AFE", "A5U")) + 
  ggsave("fa.plot_GEMIN5.pdf", width = 6, height = 5)
readRDS("encode_surf/PRPF8.results.rds") %>%
  fa.plot(plot.event = c("SE")) +
  ggsave("fa.plot_PRPF8.pdf", width = 5, height = 5)
readRDS("encode_surf/SF3B4.results.rds") %>% 
  fa.plot(plot.event = c("A3SS", "A5SS")) + 
  ggsave("fa.plot_SF3B4.pdf", width = 6, height = 5)
readRDS("encode_surf/SMNDC1.results.rds") %>% 
  fa.plot(plot.event = c("A3SS", "A5SS")) + 
  ggsave("fa.plot_SMNDC1.pdf", width = 6, height = 5)
readRDS("encode_surf/U2AF1.results.rds") %>% 
  fa.plot(plot.event = c("A3SS", "AFE")) + 
  ggsave("fa.plot_U2AF1.pdf", width = 6, height = 5)
readRDS("encode_surf/U2AF2.results.rds") %>% 
  fa.plot(plot.event = c("A3SS", "AFE")) + 
  ggsave("fa.plot_U2AF2.pdf", width = 6, height = 5)
readRDS("encode_surf/UCHL5.results.rds") %>% 
  fa.plot(plot.event = c("IAP", "TAP")) + 
  ggsave("fa.plot_UCHL5.pdf", width = 6, height = 5)


## ------ _ # of event types / RBP ------
## # detected of RBPs
fit %>% 
  filter(padj < 0.05) %>% 
  dplyr::select(factor, event) %>% 
  unique() %>%
  dplyr::count(factor)

## hist: # of associated ATR event types
g1 <- fit %>% 
  mutate(associated = padj < 0.05) %>%
  left_join(rownames_to_column(functions, var = "factor"), by = "factor") %>%
  group_by(factor, Function) %>% 
  summarise(n = length(unique(event[associated]))) %>%
  ggplot(aes(n, fill = factor(n))) + 
  geom_bar(alpha = .9, color = "grey50", size = 0.2, show.legend = F) +
  scale_fill_brewer(palette = "BrBG", direction = 1) + 
  labs(x = "# of associated event types", y = "# of RBPs") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  ggsave("bar_cnt_functions_RBP.pdf", width = 3.5, height = 3)


## #RBP at each event-feature
fit %>% 
  dplyr::filter(padj < alpha) %>% 
  group_by(factor, event, feature) %>%
  summarise(wt = 1) %>% 
  reshape2::acast(event ~ feature, sum, value.var = "wt")
#      up3 up2 up1 bd1 bd2 dn1 dn2 dn3
# SE    11  21  18  14  12  16  20  15
# RI     0   0  17  23  17  22   0   0
# A3SS  14  23  18  12  17  18   0   0
# A5SS   0   0  11  14  12  18  18  16
# AFE    0   0   6  10  18  20  22  18
# A5U    0   0  26  18  17  21   0   0
# IAP   10  18  15  14  13   9   0   0
# TAP    0   0  21  21  22  27   0   0


## #RBP each event type
fit %>% 
  dplyr::filter(padj < alpha) %>% 
  group_by(event) %>%
  summarise(n = length(unique(factor))) 
#   event     n
# 1 SE       26
# 2 RI       26
# 3 A3SS     29
# 4 A5SS     24
# 5 AFE      27
# 6 A5U      34
# 7 IAP      23
# 8 TAP      35

## #RBPs bind distal features
fit %>% filter(padj < alpha, feature %in% c("up3", "up2", "dn2", "dn3")) %>%
  summarize(distal = length(unique(factor)))
# 31

## ------ _ feature abundance ------

## chisq test of abundance
## (i) I(a location feature is associated with ATR) 
## (ii) I(a location feature is the tested feature x), x = alpha, beta, ...

abundance.test <- lapply(surf.events, function(e) {
  lapply(surf.features, function(f) {
    fit1 <- filter(fit, event == e) %>% 
      # keep only associated RBP (?)
      group_by(factor) %>%
      filter(any(padj < .05)) %>%
      ungroup() %>%
      mutate(associated = padj < .05, 
             tested = feature == f)
    if (all(!fit1$tested)) return(NULL)
    cst <- chisq.test(x = fit1$associated, y = fit1$tested)
    data.frame(event = factor(e, surf.events), 
               feature = factor(f, surf.features),
               cor = cor(x = fit1$associated, y = fit1$tested), 
               statistic = cst$statistic, 
               p.value = cst$p.value)
  }) %>% bind_rows()
}) %>% bind_rows()

pval <- reshape2::acast(abundance.test, event ~ feature, value.var = "p.value")
corr <- reshape2::acast(abundance.test, event ~ feature, value.var = "cor")
ifelse(pval > 0.05 | is.na(pval) | corr < 0, "", "*")
filter(abundance.test, event == "SE", feature == "up2") ## 0.0004881359
filter(abundance.test, event == "A3SS", feature == "up2") ## 0.0001788105


## headmap: normalized abundance
g2 <- fit %>% 
  mutate(logp = -log10(padj)) %>% 
  filter(padj < .05) %>%
  group_by(factor, event) %>% 
  mutate(wt = logp/sum(logp)) %>% 
  group_by(event, feature) %>% 
  summarise(wt = sum(wt)) %>% 
  group_by(event) %>% 
  mutate(scaleFactor = mean(wt) * 5.5, ## mean(group_size) = 5.5
         wt = wt / scaleFactor, 
         feature.type = as.integer(feature) %>% 
           cut(breaks = c(0, 3, 5, Inf), 
               labels = c("upstream", "body", "downstream"))) %>% 
  left_join(abundance.test, by = c("event", "feature")) %>%
  ggplot(aes(feature, event, fill = wt)) +
  geom_tile(color = "grey50") + 
  geom_text(aes(label = ifelse(p.value < 0.05 & cor > 0, "*", "")), 
            color = "#ffffd9", size = 7, alpha = .9, nudge_y = -0.2) + 
  labs(x = "location feature") + 
  scale_x_discrete(breaks = surf.features, labels = greek.features) +
  scale_y_discrete(limits = rev(surf.events)) +
  scale_fill_distiller(name = "abundance", palette = "GnBu", direction = 1, 
                       guide = guide_colorbar(frame.colour = "grey60")) + 
  facet_grid(cols = vars(feature.type), scales = "free_x", space = "free_x") + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_blank()) + 
  ggsave("heatmap_feature_dist.pdf", width = 4.5, height = 3)
g0 <- get_legend(g2) 
g2 <- g2 + theme(legend.position = "none")

G <- align_plots(g1, g2, align = "hv", axis = "lbrt") 
plot_grid(G[[1]], G[[2]], g0, labels = c("a", "b"), 
          ncol = 3, rel_widths = c(2, 2, .7)) + 
  ggsave("combo_fat_aggregate.pdf", width = 7, height = 3)


## ------ _ bubble plot ------ 
# min.percent = 60
# high.functional = 2.5
fit %>% 
  mutate(logp = -log10(padj)) %>%
  filter(padj < alpha) %>%
  group_by(factor, feature, event) %>%
  summarise(functional = abs.max(logp * ifelse(functional == "inclusion", 1, -1)),
            functional = sign(functional) * (abs(functional) ^ (1/2.5)),
            percent = max(logp) / sum(logp) * 100) %>% 
  ungroup() %>%
  mutate(factor = factor(factor, rev(levels(fit$factor))),
         group = functions[as.character(factor), "Function"], 
         group = replace(group, grep("&", group), "Other")) %>% 
  ggplot(aes(feature, factor, fill = functional, size = percent)) + 
  geom_point(colour = "grey70", pch = 21) +
  geom_vline(xintercept = seq(0.5, 8.5, 1), color = "grey90") + 
  geom_hline(yintercept = seq(0.5, nlevels(fit$factor) + .5, 1), color = "grey90") + 
  scale_x_discrete(breaks = surf.features, labels = greek.features) +
  # scale_y_discrete(drop = FALSE) + ## this is for plotting 76 RBPs
  scale_fill_distiller(palette = "PiYG", 
                       guide = guide_colorbar(frame.colour = "grey70")) + 
  scale_size_continuous(range = c(.5, 2.7), breaks = c(60, 80, 100)) +
  facet_grid(rows = vars(group), cols = vars(event), 
             switch = "y", scales = "free", space = "free") +
  theme_pubr(legend = "right") + 
  theme(panel.spacing.y = unit(.1, "lines"),
        strip.placement = "outside", 
        axis.text = element_text(size = 9), 
        axis.title = element_blank()) + 
  ggsave("heatmap_54rbp.pdf", width = 8, height = 7)


## ------ _ FDR estimate ------- 
## overall FDR by scrambling DR labels
fit.scrambled <- readRDS("encode_surf/all.scrambled.results.rds")
lapply(c(0.01, 0.05, 0.1), function(alpha) {
  group_by(fit.scrambled, time) %>% 
    summarise(alpha = alpha, FP = sum(padj < alpha), n = n())
}) %>% 
  list_rbind %>% 
  group_by(alpha) %>% 
  summarize(mean = mean(FP), sd = sd(FP), n = max(n))
#   alpha     mean       sd
#   <dbl>    <dbl>    <dbl>
# 1  0.01 0.000223 0.000270
# 2  0.05 0.00118  0.000727
# 3  0.1  0.00260  0.00116

## FDR stratefied by event
lapply(c(0.01, 0.05, 0.1), function(alpha) {
  group_by(fit.scrambled, time, event) %>% 
    summarise(alpha = alpha, FDR = mean(padj < alpha))
}) %>% 
  list_rbind %>% 
  group_by(alpha, event) %>% 
  summarize(mean = mean(FDR), sd = sd(FDR)) %>% 
  ggplot(aes(alpha, mean, color = event)) +
  geom_point(alpha = 0.9) + 
  geom_line(alpha = 0.8) + 
  scale_color_manual(values = setNames(colrs, events)) + 
  theme_pubr(legend = "right") + 
  ggsave("fdr_scrambled_event.pdf", width = 4, height = 5)


## ------ _ inferred features ------ 
tgWin <- readRDS("encode_surf/inferred_feature.rds")

## output GTF file 
export(tgWin, "encode_surf/inferred_feature.gtf")


## how many features are used by more than 2 RBPs
mcols(tgWin) %>% as.data.frame %>% 
  group_by(event_id, feature) %>% 
  summarise(n = length(unique(factor))) %>% 
  ungroup() %>%
  summarise("n>1" = sum(n > 1), len = n())
#   `n>1`   len
# 1 25006 71667

## ------ __ # feature / event type ------
anno_event <- readRDS("encode_surf/GRCh38.v24.filtered.event.rds")
n.feature <- setNames(c(8,4,6,6,6,4,6,4), 
                      c("SE","RI","A3SS","A5SS","AFE","A5U","IAP","TAP"))
dat0 <- table(type = anno_event$event_name) %>% data.frame %>% 
  mutate(N = Freq * n.feature[as.character(type)]) %>% 
  dplyr::select(-Freq)
dat <- mcols(tgWin) %>% as.data.frame() %>%
  group_by(event_id, feature, event_name) %>% 
  summarise(n = n()) %>% 
  group_by(event_name) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(event_name = factor(event_name, surf.events)) %>% 
  left_join(dat0, by = c("event_name" = "type")) %>%
  mutate(percent = scales::percent(n / N))
#   event_name     n      N percent
# 1 SE          9516 172792 5.51%
# 2 RI          7238  48832 14.82%
# 3 A3SS        8850 114456 7.73%
# 4 A5SS        6756  98682 6.85%
# 5 AFE         8075 192828 4.19%
# 6 A5U        11106 139884 7.94%
# 7 IAP         6372 178488 3.57%
# 8 TAP        13754 134760 10.21%

## ATR event dist: events covered inferred feature 
ggplot(dat, aes(event_name, n, fill = event_name, label = percent)) + 
  geom_bar(stat = "identity", width = .7, color = "grey60", size = .2, show.legend = F) + 
  labs(y = "# of SURF-inferred features") +
  geom_text(vjust=-.2, size = 2.5, color = "grey40") + 
  scale_fill_manual(values = c(setNames(surf.colors, surf.events), 
                               "Not inferred" = "grey60")) + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) +
  ggsave("bar_inferred_feature_dist.pdf", width = 4, height = 2.5)


## width
tgWin %>% reduce %>% width %>% sum # 5536468 

## ------ __ # RBPs / event type ------
## set data
featureSetData <- mcols(tgWin) %>% as.data.frame %>% 
  dplyr::count(factor, event_name) %>% 
  dplyr::rename("event" = "event_name") %>%
  mutate(set = paste0(factor, ".", event)) %>% 
  dplyr::select(set, everything()) %>% 
  group_by(event) %>% 
  mutate(label = paste0(event, " (", n(), ")")) %>% 
  ungroup()

label <- mcols(tgWin) %>% as.data.frame() %>%
  mutate(event = factor(event_name, surf.events)) %>%
  group_by(event) %>% 
  summarise(count = length(unique(factor))) %>% 
  mutate(label = paste0(event, " (", count, ")"))
label 
#   event count label
# 1 SE       26 SE (26)
# 2 RI       26 RI (26)
# 3 A3SS     29 A3SS (29)
# 4 A5SS     24 A5SS (24)
# 5 AFE      27 AFE (27)
# 6 A5U      32 A5U (32)
# 7 IAP      23 IAP (23)
# 8 TAP      35 TAP (35)

## ------ __ total width / event type ------
tgWin %>% split(tgWin$event_name) %>% 
  reduce %>% width %>% lapply(sum) %>% unlist
#      SE      RI    A3SS    A5SS     AFE     A5U     IAP     TAP
# 1131471  565069  942151  683510 1021678  833232  777344 1170293

## density: total width per set
g1 <- tgWin %>% split(paste(.$factor, .$event_name)) %>% 
  reduce %>% width %>% lapply(sum) %>% unlist %>% 
  data.frame(factor = sub(" .*", "", names(.)), 
             event = sub(".* ", "", names(.)) %>% factor(surf.events), 
             width = .) %>% 
  ggplot(aes(width)) + 
  geom_line(stat = "density", alpha = .7) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Total width of inferred features", 
       color = "event (# of RBPs)") + 
  theme_bw()

## box: total width per RBP, stratefied by event type
g2 <- tgWin %>% 
  split(paste0(.$factor, ".", .$event_name)) %>% 
  reduce %>% width %>% lapply(sum) %>% unlist() %>%
  data.frame(width = .) %>% 
  rownames_to_column("set") %>%
  left_join(txSetData, by = "set") %>%
  ggplot(aes(event, width, fill = label)) + 
  geom_boxplot(width = .7, alpha = .9, color = "grey40", size = 0.3) +
  scale_y_continuous(trans = "log10") + 
  scale_fill_manual(values = setNames(surf.colors, levels(txSetData$label))) +
  labs(#x = "ATR event", 
       y = "Total width of inferred features", 
       fill = "ATR event\n(# of RBPs)") + 
  coord_flip() +
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggarrange(g1, g2, ncol = 2, align = "hv", labels = "auto",
          common.legend = T, legend = "right") +
  ggsave("combo_inferred_feature_width_dist.pdf", width = 8.6, height = 3.5)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ motif analysis ------- 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------- _ extract sequences -------
library(BSgenome.Hsapiens.UCSC.hg38)

## feature set split
targetFeature = readRDS("io/encode_surf/inferred_feature.rds")
split1 <- targetFeature %>% split(.$factor) ## 52
split2 <- targetFeature %>% split(.$event_name)
split3 <- targetFeature %>%
  split(paste0(.$event_name, ".", .$feature)) ## 44
split4 <- targetFeature %>%
  split(paste0(.$factor, ".", .$event_name)) ## 222
split5 <- targetFeature %>%
  split(paste0(.$factor, ".", .$event_name,".", .$functional)) ## 365
grl <- c(split1, split2, split3, split4, split5)

#   52(RBP)
# + 8(ATR)
# + 44(ATR x feature)
# + 222(RBP x ART)
# + 365(RBP x ART x functional)
# = 691

## get genomic sequences
for (i in seq_along(grl)) {
  getSeq(BSgenome.Hsapiens.UCSC.hg38, grl[[i]]) %>%
    RNAStringSet() %>%
    writeXStringSet(paste0("io/encode_surf/",
                           names(grl)[i],
                           ".feature.fasta"))
}

grl <- c(split1, split2, split3, split4)
write.table(names(grl), "featureSet.txt",
            quote = F, row.names = F, col.names = F)

nm <- sub(".inclusion|.exclusion","",names(split5))
which(table(nm) == 2) %>% names %>%
  write.table("featureSet2.txt", quote = F,
              row.names = F, col.names = F)
## 143 RBPxATR combinations

## overlapping
nm <- sub(".inclusion|.exclusion","",names(split5))
ss <- which(table(nm) == 2) %>% names()
percOverlaps <- sapply(ss, function(s) {
  x = split5[[paste0(s, ".inclusion")]]
  y = split5[[paste0(s, ".exclusion")]]
  sum(!!countOverlaps(x, y)) / (length(x) + length(y))
}) %>% setNames(ss)
summary(percOverlaps)

targetFeature = readRDS("io/encode_surf/inferred_feature.rds")
nm <- paste0(targetFeature$factor, ".", targetFeature$event_name) %>% unique
info <- data.frame(
  name = nm, 
  factor = sub("\\..*", "", nm), 
  event = factor(sub(".*\\.", "", nm), surf.events), 
  row.names = nm, 
  stringsAsFactors = F
)

## load motifs: RBP x ATR
motifs <- list() 
err <- c("HNRNPM.TAP", "PUM2.TAP", "EIF3G.A5SS") ## these three set has only one inferred feature
for (s in setdiff(nm, err)) {
  filename <- paste0("io/encode_surf/", s, ".meme/meme.txt")
  motif <- read_meme(filename) %>% 
    filter_motifs(eval = 0.01)
  if (length(motif)) {
    motifs[[s]] <- motif[[1]]
    motifs[[s]]@name = s
  }
}


## how many RBP has identified motif?
info[names(motifs),] %>% dplyr::count(factor) 
## 30 RBPs

## how about event type?
info[names(motifs),] %>% dplyr::count(event) 

g1 <- view_motifs(motifs$HNRNPC.TAP) +
  scale_y_continuous(limits = c(0,2)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank()) +
  ggsave("motif_HNRNPC.TAP.pdf", width = 3, height = 1.5)
g2 <- view_motifs(motifs$HNRNPC.SE) +
  scale_y_continuous(limits = c(0,2)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank()) +
  ggsave("motif_HNRNPC.SE.pdf", width = 3, height = 1.5)
plot_grid(g1, g2, nrow = 1, 
          align = "hv", axis = "rlbt") + 
  ggsave("motif_HNRNPC.pdf", width = 6, height = 1.5)

view_motifs(motifs[c("TIA1.AFE","TIA1.IAP","TIA1.TAP")]) + 
  theme_classic() +
  theme(axis.title.x = element_blank()) + 
  ggsave("motif_TIA1.pdf", width = 5, height = 5)

## ------ _ diversity ------- 
sim.pwm <- compare_motifs(motifs, method = "KL") %>%
  reshape2::melt(varnames = c("m1", "m2"), 
                 value.name = "sim") %>% 
  left_join(info, c("m1" = "name")) %>% 
  left_join(info, c("m2" = "name")) %>% 
  filter(m1 != m2, factor.x == factor.y) %>% 
  group_by(factor.x) %>% 
  summarise(n = length(unique(event.x)), 
            min = min(sim), 
            max = max(sim), 
            mean = mean(sim)) %>% 
  ungroup() %>%
  arrange(mean) %>% 
  mutate(label = paste0(factor.x, " (", n, ")"), 
         label = factor(label, label)) 

sim.pwm %>% 
  ggplot(aes(label, mean, color = summary)) + 
  geom_pointrange(aes(ymin = min, ymax = max), 
                  color = "#E69F00", alpha = .9) + 
  labs(x = "RBP (# of ATR-specific motifs)", 
       y = "Kullback-Leibler divergence") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  ggsave("scatter_sim_pwm.pdf", width = 5, height = 3.4)

## ------ _ clustering ------- 
## dendro: PWM similarity
sim.pwm <- compare_motifs(motifs, method = "KL")
tree.sim <- hclust(as.dist(sim.pwm), method = "ward.D2") %>% 
  ape::as.phylo() %>%
  ggtree(layout = "fan",
         open.angle = 90) %>% 
  rotate_tree(-45)
tree.sim %<+% info +
  geom_tiplab2(aes(color = event)) + 
  scale_color_manual(values = surf.colors) +
  ggsave("dendro_motif_tree_KL.pdf", width = 7, height = 6)

## ------- _ gr overlaps ------- 
# how many overlaps
fs <- split(targetFeature, paste0(targetFeature$factor, ".", targetFeature$event_name))
percOverLaps <- sapply(names(motifs), function(s) {
  nm.other <- setdiff(names(motifs), s)
  event.other <- sub(".*\\.", "", nm.other)
  event.this <- sub(".*\\.", "", s)
  other <- fs[nm.other[event.other == event.this]]
  nOverlaps <- sapply(other, function(f) {
    sum(!!countOverlaps(fs[[s]], f))
  })
  mean(nOverlaps) / length(fs[[s]])
})
data.frame(name = names(percOverLaps), 
           percent = percOverLaps) %>% 
  left_join(info, "name") %>% 
  group_by(event) %>%
  mutate(label = paste0(event, " (", n(), ")")) %>% 
  ungroup() %>% 
  arrange(event) %>%
  mutate(label = factor(label, unique(label))) %>%
  ggplot(aes(event, percent, color = label)) + 
  geom_boxplot(alpha = .9, color = "grey50") + 
  geom_jitter(alpha = .8, width = 0.12) + 
  # geom_histogram(alpha = .8, color = "grey70", 
  #                size = .2, bins = 15) + 
  scale_y_continuous("average pair-wise overlaps", 
                     labels = scales::percent) + 
  scale_color_manual("event type\n(# of RBPs)", 
                     values = surf.colors) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  ggsave("hist_percent_overlap_feature.pdf", 
         width = 6, height = 2.5)

## # of overlaps location features
fs <- split(targetFeature, paste0(targetFeature$factor, ".", targetFeature$event_name))
nOverLaps <- sapply(fs[names(motifs)], function(x) {
  sapply(fs[names(motifs)], function(y) {
    sum(!!countOverlaps(x, y))
  })
})

## use proportion
# nOverLaps <- sapply(fs[names(motifs)], function(x) {
#   sapply(fs[names(motifs)], function(y) {
#     sum(!!countOverlaps(x, y)) / 
#       length(unique(c(names(x), names(y))))
#   })
# })

## cluster based on # overlaps
nInverse <- diag(sqrt(1 / diag(nOverLaps)))
dimnames(nInverse) = dimnames(nOverLaps)
sim.overlaps <- as.dist(1 - nInverse %*% nOverLaps %*% nInverse)
tree <- hclust(sim.overlaps, method = "ward.D2") %>% 
  ape::as.phylo() %>% 
  ggtree(layout = "fan", open.angle = 90, 
         branch.length = "none") 
tree %<+% info %>%
  rotate_tree(-45) +
  geom_tiplab2(aes(color = event)) + 
  scale_color_manual(breaks = surf.events,
                     values = setNames(surf.colors, surf.events)) +
  ggsave("dendro_inferredFeature_overlaps.pdf", width = 7, height = 6)

row.names <- clusterByGroup(nOverLaps, info[rownames(nOverLaps), "event"])
# col.names <- clusterByGroup(t(nOverLaps), info[colnames(nOverLaps), "event"])
dat <- nOverLaps[row.names, row.names] %>%
  reshape2::melt(varnames = c("row", "column"), 
                 value.name = "count") %>%
  left_join(info, c("row" = "name")) %>% 
  left_join(info, c("column" = "name")) %>% 
  mutate(factor.x = factor(factor.x, unique(factor.x)), 
         factor.y = factor(factor.y, rev(unique(factor.y)))) 
dat %>% 
  ggplot(aes(factor.x, factor.y, fill = count)) +
  geom_tile(color = "grey80", size = .1) +
  facet_grid(rows = vars(event.y), cols = vars(event.x), 
             scales = "free", space = "free") + 
  scale_fill_distiller(trans = "sqrt", palette = "GnBu", direction = 1,
                       breaks = c(50, 500, 1500, 3000), 
                       guide = guide_colorbar(frame.colour = "grey60")) + 
  theme_classic() +
  theme(panel.spacing = unit(.2, "lines"), 
        axis.title = element_blank(), 
        # axis.text = element_blank(),
        # axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5)) + 
  ggsave("heatmap_feature_overlaps_count.pdf", width = 11.5, height = 10.5)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ DASeq (module 3) ------- 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

targetFeature = readRDS("encode_surf/inferred_feature.rds") %>% data.frame()
study <- c("TCGA", "GTEx")
color <- c("#56B4E9", "#E69F00")
# dar <- readRDS("encode_surf/daseq.results.v1.rds")
dar <- readRDS("encode_surf/daseq.results.v2.rds")
# dar <- readRDS("encode_surf/daseq.results.gene.rds")
sampleData <- rownames_to_column(data.frame(dar@sampleData), "sample")

targetTxSets <- rowData(dar@targetAUC)$set
exprMat.tx <- readRDS("TcgaTargetGtex_rsem_isoform_tpm_laml_blood.rds")
isopct <- readRDS("TcgaTargetGtex_rsem_isopct_laml_blood.rds")
rankings.tx <- readRDS("TcgaGtex_tx_ranking_laml_blood.rds") 

targetGeneSets <- targetFeature$gene_id %>% 
  split(targetFeature[,c("factor", "event_name")]) %>%
  lapply(unique) %>% purrr::compact()
exprMat.gene <- readRDS("TcgaTargetGtex_rsem_gene_tpm_laml_blood.rds")
rankings.gene <- readRDS("TcgaGtex_gene_ranking_laml_blood.rds")

anno <- import("annotation/gencode.v24.annotation.tpm1isoPct5.gtf", "gtf")
anno.full = import("annotation/gencode.v24.annotation.gtf", "gtf")

## ------ _ gene level ------ 

## ------ __ sample stat ------ 
## # of genes in each sample
gene_protein <- anno.full$gene_id[anno.full$type == "gene" & anno.full$gene_type == "protein_coding"]
g2 <- data.frame(sample = colnames(exprMat.gene), 
                nGene = colSums(exprMat.gene[rownames(exprMat.gene) %in% gene_protein,] > log2(0.001))) %>% 
  left_join(sampleData, by = "sample") %>% 
  ggplot(aes(nGene, y = ..density.., fill = condition)) + 
  geom_histogram(color = "grey70", size = .2, bins = 30,
                 alpha = .8, position = "identity") + 
  labs(x = "# of protein-coding genes", fill = "study") + 
  scale_fill_manual(values = setNames(color, study)) + 
  theme_bw() + 
  ggsave("hist_n_genes_tcga_gtex.pdf", width = 4, height = 2.8)

## ------ __ set stat ------ 
gene_protein <- anno.full$gene_id[anno.full$type == "gene" & anno.full$gene_type == "protein_coding"]
geneSetData <- lapply(targetGeneSets, function(s) {
  data.frame(size = length(s), 
             perc_protein = mean(s %in% gene_protein))
}) %>% 
  bind_rows(.id = "set") %>% 
  mutate(factor = sub("\\..*", "", set),
         event = factor(sub(".*\\.", "", set), surf.events)) %>% 
  arrange(event, factor) %>%
  group_by(event) %>% 
  mutate(label = paste0(event, " (", n(), ")")) %>%
  ungroup() %>% 
  mutate(label = factor(label, unique(label)))

## how many genes
summary(geneSetData$size) 
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     1.0    45.0   152.5   271.6   397.8  1556.0

geneSetData$size %>% 
  cut(c(0, 50, 100, 150, Inf)) %>% 
  table(size = .) 
#   (0,50]  (50,100] (100,150] (150,Inf]
#       57        34        19       112

## how many are protein coding 
summary(geneSetData$perc_protein) ## mean=97.34%
arrange(geneSetData, perc_protein)

## count target genes grouped by ATR event
targetGeneSetsByEvent <- targetFeature$gene_id %>%
  split(targetFeature[, "event_name"]) %>%
  lapply(unique) %>% purrr::compact()
cntGene <- sapply(targetGeneSetsByEvent, function(x) {
  sapply(targetGeneSetsByEvent, function(y){
    length(intersect(x, y))
  })}) %>% 
  reshape2::melt(value.name = "count") %>% 
  mutate(Var1 = factor(Var1, surf.events), 
         Var2 = factor(Var2, rev(surf.events)))


## ------ __ AUC ------ 

## heat map raw auc
mat <- calculateAUC(targetGeneSets, rankings.gene) %>% getAUC()
clustered_targetSet = sub(".*\\.", "", rownames(mat)) %>% 
  factor(surf.events) %>% clusterByGroup(x = mat)
clustered_sample <- t(mat) %>% 
  clusterByGroup(sampleData[colnames(mat), "condition"])
g <- reshape2::melt(mat, value.name = "AUC") %>% 
  left_join(geneSetData, by = "set") %>%
  mutate(condition = sampleData[sample, "condition"]) %>%
  ggplot(aes(sample, set, fill = AUC)) +
  geom_raster() + 
  scale_x_discrete(breaks = clustered_sample) +
  scale_y_discrete(breaks = clustered_targetSet) +
  scale_fill_distiller(name = "AUC", palette = "GnBu", direction = 1,
                       guide = guide_colorbar(frame.colour = "grey60")) + 
  facet_grid(rows = vars(event), cols = vars(condition), 
             scales = "free", space = "free") + 
  labs(x = "Sample", y = "Target gene set") +
  theme_bw() +
  theme(panel.spacing = unit(.2, "lines"), 
        panel.border = element_rect(color = "grey70"),
        panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank()) + 
  ggsave("heatmap_auc_geneset.pdf", width = 7, height = 7)

## target gene sets AUC (tidy)
auc.target <- calculateAUC(targetGeneSets, rankings.gene) %>% 
  tidyAUC() %>% mutate(type = "target")

## randomly sampled gene set, matching sizes and ranks
gene_id <- assays(rankings.gene)$ranking %>% 
  rowMeans() %>% sort(decreasing = F) %>% names
randomGeneSets <- lapply(targetGeneSets, function(s) {
  ranks <- match(s, gene_id) 
  stepSize = 1e3L
  binCount <- table(ranks %/% stepSize)
  naCount <- sum(is.na(ranks))
  lapply(seq_along(binCount), function(i) {
    gene_id[(stepSize * (i - 1) + 1):(stepSize * i)] %>% 
      sample(binCount[i])
  }) %>% unlist() %>% 
    append(sample(gene_id, naCount))
  # r <- 1.4 * mean(match(s, gene_id), na.rm = T)
  # sample(gene_id[1:r], length(s))
})
auc.random <- calculateAUC(randomGeneSets, rankings.gene) %>% 
  tidyAUC() %>% mutate(type = "random")

## compare AUC, TCGA vs GTEx
dat <- bind_rows(auc.target, auc.random) %>% 
  left_join(sampleData, by = "sample") %>% 
  group_by(set, condition, type) %>% 
  summarise(AUC = mean(AUC, na.rm = T)) %>% 
  dplyr::ungroup() %>%
  reshape2::dcast(set + type ~ condition, value.var = "AUC") %>% 
  mutate(base = (TCGA + GTEx) / 2, 
         diff = TCGA - GTEx) %>% 
  left_join(geneSetData, by = "set")

## scatter: TCGA vs GTEx
## We have tried coloring dots with event, but there is little to see
g1 <- ggplot(dat, aes(GTEx, TCGA, color = type)) + 
  geom_point(alpha = .3) + 
  geom_abline(slope = 1, linetype = 2, intercept = 0, alpha = .6) + 
  geom_rug(size = .2, alpha = .8) + 
  labs(x = "activity level in GTEx", 
       y = "activity level in TCGA", 
       color = "gene set") + 
  scale_color_manual(values = c(target = "#4eb3d3", 
                                random = "grey70")) +
  theme_bw() + 
  ggsave("scatter_geneset_auc.pdf", width = 4, height = 2.8)
g1 <- g1 + theme(legend.position = "none")

## box: diff
dat %>% 
  filter(type == "target") %>% 
  ggplot(aes(event, diff, fill = event, label = label)) + 
  geom_boxplot(alpha = .9, color = "grey50", varwidth = T) + 
  geom_hline(yintercept = 0, color = "grey40", linetype = 2) +
  scale_fill_manual(labels = levels(geneSetData$label), 
                    values = setNames(surf.colors, surf.events)) + 
  labs(y = bquote(Delta*"AUC (TCGA vs GTEx)"), 
       fill = "ATR event\n(# of gene sets)") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank()) + 
  ggsave("box_geneset_auc.pdf", width = 5, height = 3)


## compare dispersions: TCGA vs GTEx
auc.target %>% 
  left_join(sampleData, by = "sample") %>% 
  group_by(set, condition) %>% 
  summarise(sd = sd(AUC, na.rm = T)) %>% 
  dplyr::ungroup() %>%
  reshape2::dcast(set ~ condition, value.var = "sd") %>% 
  left_join(geneSetData, by = "set") %>%
  ggplot(aes(TCGA, GTEx, color = event)) +
  geom_point(alpha = .5, show.legend = F) +
  geom_abline(slope = 1, color = "grey60", linetype = 2) +
  scale_x_continuous(limits = c(0,.1)) + 
  scale_y_continuous(limits = c(0,.1)) + 
  labs(x = "deviation of AUC (TCGA)", y = "deviation of AUC (GTEx)") +
  scale_color_manual(labels = levels(geneSetData$label), 
                     values = setNames(surf.colors, surf.events)) + 
  facet_wrap(~ label, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  ggsave("scatter_geneset_auc_sd.pdf", width = 7.2, height = 4.2)


## ------ _ transcript level ------ 

## ------ __ set stat ------ 
## set Data
txSetData <- lapply(targetTxSets, function(s) {
  data.frame(size = length(s))
}) %>% 
  bind_rows(.id = "set") %>% 
  mutate(factor = sub("\\..*", "", set),
         event = factor(sub(".*\\.", "", set), surf.events)) %>% 
  arrange(event, factor) %>%
  group_by(event) %>% 
  mutate(label = paste0(event, " (", n(), ")")) %>%
  ungroup() %>% 
  mutate(label = factor(label, unique(label)))

## size of target sets
txSetData$size %>% cut(c(0, 50, 100, 150, Inf)) %>% table(size = .) 
#  (0,50]  (50,100] (100,150] (150,Inf]
#      57        33        20       112
summarise(txSetData, total = n(), 
          ">50" = sum(size > 50), 
          ">100" = sum(size > 100))
# total   >50   >100 
#   222   165    132


## how many element are unique to each set? -- very few
counter <- unlist(targetTxSets) %>% table 
uniques <- sapply(targetTxSets, function(x) {
  return(sum(counter[x] == 1))
})
hist(uniques); dev.off()

## how much do set category share transcripts
targetTxSetsByEvent <- targetFeature$transcript_id %>%
  split(targetFeature[, "event_name"]) %>%
  lapply(unique) %>% purrr::compact()
cntTx <- sapply(targetTxSetsByEvent, function(x) {
  sapply(targetTxSetsByEvent, function(y){
    length(intersect(x, y))
  })}) %>% 
  reshape2::melt(value.name = "count") %>% 
  mutate(Var1 = factor(Var1, surf.events), 
         Var2 = factor(Var2, rev(surf.events))) 

## shared target gene/tx by event types
g1 <- bind_rows(gene = cntGene, 
          transcript = cntTx, 
          .id = "target") %>%
  ggplot(aes(Var1, Var2, fill = count)) + 
  geom_tile(color = "grey50") + 
  scale_fill_distiller(name = "# of targets", 
                       palette = "GnBu", direction = 1, 
                       guide = guide_colorbar(frame.colour = "grey60")) + 
  facet_wrap(~ target, nrow = 1) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title = element_blank()) + 
  ggsave("heatmap_shared_target_by_event.pdf", 
         width = 7, height = 3.4)

## box: target set sizes
g2 <- bind_rows("gene set" = geneSetData, 
                "transcript set" = txSetData, 
                .id = "target") %>% 
  ggplot(aes(event, size, fill = label)) + 
  geom_boxplot(alpha = .9, varwidth = T, color = "grey50", 
               outlier.colour = "grey70") + 
  labs(y = "# of targets", fill = "event type\n(# of sets)") +
  scale_y_continuous(trans = "sqrt", breaks = c(10, 100, 500, 1500)) + 
  scale_fill_manual(values = surf.colors) + 
  facet_wrap(~ target, nrow = 1) + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  ggsave("box_size_target_set_by_event.pdf", width = 7, height = 3)

plot_grid(g1, g2, align = "hv", axis = "rlbt", ncol = 1, 
          labels = "auto", rel_heights = c(1,1)) + 
  ggsave("combo_target_count.pdf", width = 8, height = 7)


## ------ __ AUC (not used) ------ 

## target gene sets AUC
auc.target <- calculateAUC(targetTxSets, rankings.tx) %>% tidyAUC() 

## randomly sampled gene set, matching sizes and ranks
tx_id <- rowMeans(assays(rankings.tx)$ranking) %>% 
  sort(decreasing = F) %>% names
randomTxSets <- lapply(targetTxSets, function(s) {
  r <- 1.4 * mean(match(s, tx_id), na.rm = T)
  sample(tx_id[1:r], length(s))
})
auc.random <- calculateAUC(randomTxSets, rankings.tx) %>% tidyAUC()

## base line
# mead(auc.random$AUC) ## 0.1133494

## contrast auc
dat <- bind_rows("target" = auc.target, 
                 "random" = auc.random, 
                 .id = "type") %>% 
  left_join(sampleData, by = "sample") %>% 
  group_by(set, condition, type) %>% 
  summarise(AUC = mean(AUC, na.rm = T)) %>% 
  dplyr::ungroup() %>%
  reshape2::dcast(set + type ~ condition, value.var = "AUC") %>% 
  mutate(diff = TCGA - GTEx) %>% 
  left_join(txSetData, by = "set")

## scatter: TCGA vs GTEx
dat %>%
  ggplot(aes(GTEx, TCGA, color = type)) + 
  geom_point(alpha = .3) + 
  geom_abline(slope = 1, linetype = 2, intercept = 0, alpha = .6) + 
  labs(x = "activity level in GTEx", y = "activity level in TCGA", color = "set") + 
  scale_color_manual(values = c(target = "#56B4E9", random = "grey70")) +
  theme_bw() + 
  ggsave("scatter_txset_auc.pdf", width = 4, height = 3)

## ------ __ isoPct ------ 
getIsoPct <- function(tx_id) {
  dat.isopct <- isopct[tx_id,] %>% 
    rownames_to_column("transcript_id") %>%
    pivot_longer(-transcript_id, "sample", values_to = "isopct") %>%
    left_join(sampleData, by = "sample") %>%
    group_by(transcript_id, condition) %>% 
    summarise(isopct = mean(isopct, na.rm = T)) %>% 
    ungroup() %>%
    pivot_wider(names_from = condition, 
                names_prefix = "isopct.",
                values_from = isopct) %>% 
    mutate(isopct.diff = isopct.TCGA - isopct.GTEx) 
  dat.tpm <- exprMat.tx[tx_id,] %>%
    rownames_to_column("transcript_id") %>%
    pivot_longer(-transcript_id, "sample", values_to = "tpm") %>%
    left_join(sampleData, by = "sample") %>%
    group_by(transcript_id, condition) %>% 
    summarise(tpm = mean(tpm, na.rm = T)) %>% 
    ungroup() %>%
    pivot_wider(names_from = condition, 
                names_prefix = "tpm.", 
                values_from = tpm) %>%
    mutate(tpm.diff = tpm.TCGA - tpm.GTEx)
  
  left_join(dat.isopct, dat.tpm, by = "transcript_id")
}

summariseIsoPct <- function(tbl) {
  if (nrow(tbl) <= 4) t = NULL else 
    t <- t.test(tbl$isopct.TCGA, 
                tbl$isopct.GTEx, 
                paired = TRUE)
  ## mean isoPct within condition
  tbl %>% summarise(
    TCGA = mean(isopct.TCGA), 
    GTEx = mean(isopct.GTEx), 
    diff = mean(isopct.TCGA - isopct.GTEx), 
    base = mean((isopct.TCGA + isopct.GTEx) / 2)
  ) %>% mutate(
    df = t$parameter, 
    stat = t$stat, 
    p.value = t$p.value
  )
}


## target set isopct
isopct.target <- targetTxSets %>% 
  lapply(getIsoPct) %>% 
  lapply(summariseIsoPct) %>% 
  bind_rows(.id = "set") %>% 
  mutate(padj = p.adjust(p.value, "fdr")) %>% 
  left_join(txSetData, by = "set") 


## # of significance by paired t-test
sum(isopct.target$padj < .05, na.rm = T) ## 93
sum(isopct.target$padj < .05 & isopct.target$df > 100, na.rm = T) ## 75

## randomly sampled tx set
# tx_id <- rowMeans(assays(rankings.tx)$ranking) %>% 
#   sort(decreasing = F) %>% names
isopct.random <- lapply(seq_len(5), function(i) {
  mapply(function(ts, gs) {
    anno$transcript_id[anno$gene_id %in% gs] %>% 
      unique() %>% sample(length(ts))
  }, targetTxSets, targetGeneSets, SIMPLIFY = F) %>%
    lapply(getIsoPct) %>% 
    lapply(summariseIsoPct) %>% 
    bind_rows(.id = "set") %>% 
    mutate(padj = p.adjust(p.value, "fdr")) %>% 
    left_join(txSetData, by = "set") %>%
    mutate(control = i)
}) %>% bind_rows()

## scatter: base vs diff
isopct.target %>% 
  ggplot(aes(base, diff, color = event)) + 
  geom_point(alpha = .3) + 
  geom_hline(yintercept = mean(isopct.random$diff), 
             linetype = 2, alpha = .6) + 
  labs(x = "mean IsoPct", 
       y = bquote(Delta*"IsoPct (TCGA vs GTEx)"), 
       color = "transcript set") + 
  scale_y_continuous(limits = c(-10, 10)) + 
  scale_color_manual(values = setNames(surf.colors, surf.events)) +
  theme_bw() + 
  ggsave("scatter_txset_isopct.pdf", width = 4, height = 2.8)

## box: diff
g2 <- bind_rows("target" = isopct.target, 
                "random" = isopct.random,
                .id = "type") %>% 
  mutate(group = paste(event, control)) %>%
  filter(control == 1 | type == "target") %>% ## this is for main figure
  ggplot(aes(event, diff, fill = type, 
             group = paste(event, control))) + 
  geom_boxplot(alpha = .9, varwidth = T, color = "grey50", 
               outlier.colour = "grey70") + 
  geom_hline(yintercept = 0, color = "grey40", linetype = 2) +
  scale_y_continuous(limits = c(-10, 10)) + 
  scale_fill_manual(values = c(target = "#4eb3d3", random = "grey70")) +
  labs(y = bquote(Delta*"IsoPct (TCGA vs GTEx)"), 
       fill = "transcript set") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) +
  ggsave("box_isoPct_diff.pdf", width = 8, height = 3.5)
g0 = get_legend(g2)
g2 = g2 + theme(legend.position = "none")

G <- align_plots(g1, g2, align = "hv", axis = "rlbt") 
plot_grid(G[[1]], G[[2]], g0, labels = c("a", "b"), 
          nrow = 1, rel_widths = c(2.7,3.6,1)) + 
  ggsave("combo_target_gene_tx.pdf", width = 7.5, height = 2.7)

## ------ ___ ma plot ------ 
## scatter: DeltaIsoPct of transcripts
selectedSet = "SRSF1.IAP"
selectedSet = c("PUM1.A3SS", "HNRNPC.A5SS", 
                "FXR1.RI", "IGF2BP2.TAP")
g1 <- split(selectedSet, selectedSet) %>% 
  lapply(function(s) {
    anno$gene_id %in% targetGeneSets[[s]] %>% 
      anno$transcript_id[.] %>% 
      unique() %>% getIsoPct() %>% 
      mutate(transcript = transcript_id %in% targetTxSets[[s]], 
             transcript = ifelse(transcript, "target", "random")) %>% 
      mutate(base = (tpm.GTEx + tpm.TCGA) / 2,
             change = isopct.TCGA - isopct.GTEx) 
  }) %>% bind_rows(.id = "set") %>% 
  left_join(txSetData, by = "set") %>% 
  mutate(label = paste0(set, " (", size, ")")) %>% 
  ggplot(aes(base, change)) +
  geom_point(aes(shape = transcript), color = "#4eb3d3", alpha = .2) +
  stat_summary(aes(x = 0, yintercept = ..y.., 
                   color = transcript, 
                   linetype = transcript), 
               fun.y = mean, geom = "hline") +  
  geom_rug(size = .1, alpha = .5, color = "grey70", show.legend = F) + 
  labs(x = bquote("base expression (log"[2]~"TPM)"), 
       y = bquote(Delta*"IsoPct (TCGA vs. GTEx)")) + 
  scale_color_manual(values = c("target" = "#0868ac", "random" = "grey42")) +
  scale_shape_manual(values = c("target" = 19, "random" = NA)) + 
  facet_wrap(~ label, ncol = 2) +
  theme_bw() + 
  ggsave(paste0("ma_isopct_vs_tpm_tx.pdf"), width = 7, height = 6)
  # scale_y_continuous(limits = c(-20,20)) + 
  # ggsave(paste0("ma_isopct_vs_tpm_tx_", selectedSet,".pdf"),
  #        width = 4, height = 3)

## ------ ___ group by event type ------ 
targetTxSetsByEventUnique <- targetFeature$transcript_id %>%
  split(targetFeature[, "event_name"]) %>%
  lapply(unique) %>% 
  lapply(as.data.frame, stringsAsFactors = F) %>% 
  bind_rows(.id = "event") %>% 
  mutate(event = factor(event, surf.events)) %>%
  dplyr::rename("transcript_id" = "X[[i]]") %>% 
  group_by(transcript_id) %>% 
  mutate(n = n()) %>% 
  group_by(event) %>%
  mutate(label = paste0(event, " (", 
                        sum(n == 1), "/",
                        n(), ")")) %>%
  ungroup() %>%
  mutate(label = factor(label, unique(label))) %>%
  filter(n == 1) %>% 
  left_join(getIsoPct(.$transcript_id), "transcript_id") 

## pair-wise t test
pairwise.t.test(targetTxSetsByEventUnique$isopct.diff, 
                targetTxSetsByEventUnique$event, 
                p.adjust.method = "fdr") 

# Pairwise comparisons using t tests with pooled SD
# 
#        SE      RI      A3SS    A5SS    AFE     A5U     IAP
#   RI   0.9778  -       -       -       -       -       -
#   A3SS 0.0851  0.4703  -       -       -       -       -
#   A5SS 0.2419  0.6751  0.9778  -       -       -       -
#   AFE  0.0086  0.0049  6.9e-07 1.5e-05 -       -       -
#   A5U  1.0e-06 3.1e-06 2.9e-12 8.4e-10 0.2419  -       -
#   IAP  < 2e-16 < 2e-16 < 2e-16 < 2e-16 < 2e-16 < 2e-16 -
#   TAP  < 2e-16 < 2e-16 < 2e-16 < 2e-16 3.0e-10 3.1e-06 1.5e-05
# 
# P value adjustment method: holm

## box: isopct by event type
targetTxSetsByEventUnique %>% 
  ggplot(aes(label, isopct.diff, fill = label)) + 
  geom_boxplot(alpha = .8, varwidth = T, color = "grey40") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") + 
  scale_y_continuous(trans = trans_new("sqrt.new", 
                                       function(x) sign(x) * sqrt(abs(x)), 
                                       function(x) sign(x) * x^2), 
                     breaks = c(-30, -10, -2, 0, 2, 10, 30)) + 
  labs(y = bquote(Delta*"IsoPct (TCGA vs GTEx)"), 
       fill = "event type\n(# of unique /\ntotal transcripts)") +
  scale_fill_manual(values = surf.colors) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ggsave("box_isopct_unique_tx_event.pdf", width = 6, height = 3)


## ------ __ DA ------
tidy.dar <- data.frame(dar) %>% 
  mutate(baseAUC = (GTEx + TCGA) / 2,
         dAUC = GTEx - TCGA) %>% 
  left_join(txSetData, by = c("set")) 
dat1 <- filter(tidy.dar, padj < 0.1, size > 100) %>% na.omit()
dat1 %>% arrange(-abs(stat)) %>%
  dplyr::select("set", "size", "GTEx", "TCGA", 
         "stat", "p.value", "padj")

## check dist of p-value
tidy.dar %>%
  # filter(size > 100) %>%
  ggplot(aes(p.value, fill = label)) + 
  geom_histogram(alpha = .8, size = .2, bins = 11, color = "grey70") + 
  labs(x = "raw p-value", fill = "event type\n(# of sets)") + 
  scale_fill_manual(values = surf.colors) +
  theme_bw() + 
  ggsave("hist_dar_pvalue_11bins.pdf", width = 4, height = 2.8)

## ------ ___ basic count ------ 
## # of DA
sum(tidy.dar$padj < 0.05, na.rm = T) ## 19
sum(tidy.dar$padj < 0.05 & tidy.dar$size > 100, na.rm = T) ## 18
sum(tidy.dar$padj < 0.1, na.rm = T) ## 23
sum(tidy.dar$padj < 0.1 & tidy.dar$size > 100, na.rm = T) ## 19


## relevant RBP counts stratefied by event type
table(dat1$event)
#   SE   RI A3SS A5SS  AFE  A5U  IAP  TAP
#    4    1    7    7    0    0    0    0

## RBP counts stratefied by event type
table(dat1$factor)
#  AQR   FXR1 HNRNPC  PRPF8   PUM1  RBM15  SF3B4  SRSF1  U2AF2
#    2      4      2      1      1      2      1      3      3


## print DA set details
tidy.dar %>% 
  filter(size > 10) %>%
  dplyr::select(factor, event, size, GTEx, TCGA, padj) %>% 
  mutate(significance = cut(padj, 
                            c(-1, 0.001, 0.01, 0.05, 0.1, 1), 
                            c("***", "**", "*", ".", ""))) %>%
  dplyr::rename("AUC(GTEx)" = "GTEx", 
         "AUC(TCGA)" = "TCGA", 
         "adj. p-value" = "padj") %>%
  xtable::xtable(digits = 4)


## ------ ___ ma plot ------

## scatter: isopct
pickSet <- c("PUM1.A3SS", "HNRNPC.A5SS", "AQR.A5SS", "FXR1.A3SS", "FXR1.RI")
auc.control = dar@controlAUC %>% aggregateAUCbyCondition()
g2 <- tidy.dar %>% 
  filter(size > 10) %>%
  mutate(color = ifelse(set %in% dat1$set, "diff.", "equal"), 
         label = ifelse(set %in% pickSet, set, "")) %>%
  ggplot(aes(baseAUC, dAUC, color = color)) + 
  geom_point(alpha = .4) + 
  geom_hline(yintercept = mean(auc.control$diff), linetype = 2, alpha = .6) + 
  scale_y_continuous(limits = c(-0.2, 0.1)) +
  labs(x = "base activity level (mean AUC)", 
       y = bquote(Delta*"AUC (TCGA vs GTEx)"), 
       color = "activity") + 
  scale_color_manual(values = c("#fb8072", "grey60")) +
  ggrepel::geom_label_repel(
    aes(label = label), size = 2, color = "grey50", force = 6,
    label.padding = .16, box.padding = 0.4, point.padding = 0.2,
    segment.size = 0.3, segment.alpha = .8, segment.color = "grey50", 
    show.legend = F) +
  theme_bw() + 
  ggsave("ma_txset_dar.pdf", width = 3.7, height = 2.5)

plot_grid(g1, g2, nrow = 1, align = "hv", axis = "rlbt", 
          rel_widths = c(2.1, 2)) + 
  ggsave("combo_ma_SRSF1.IAP_dar.pdf", width = 7.9, height = 2.9)


