


<!-- README.md is generated from README.Rmd. Please edit that file -->

# Code for SURF Analysis of ENCODE Data

<!-- badges: start -->

<!-- [![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) -->

<!-- badges: end -->

This repository makes available the source code for our SURF paper(s).
Last updated in April 2020.

The paper presents the **S**tatistical **U**tility for **R**BP
**F**unctions (SURF) for integrative analysis of RNA-seq and CLIP-seq
data. The goal of SURF is to identify alternative splicing (AS),
alternative transcription initiation (ATI), and alternative
polyadenylation (APA) events regulated by individual RBPs and elucidate
protein-RNA interactions governing these events. We apply the SURF
pipeline to analyze 104 RBP data sets (from
[ENCODE](https://www.encodeproject.org)). Check out the browsable
results from this [shiny](http://www.statlab.wisc.edu/shiny/surf/) app\!

The current repository includes:

  - application/
      - `xena.R`: process TCGA and GTEx transcriptome data.
      - `encode_surf_one.R`: perform SURF analysis for one RBP. This is
        used for all 104 RBPs.
      - `encode_surf_summary.R`: summarize the SURF results, including
        all the statistics and plots reported in the paper.
  - simulation/
      - `other_simulation.sh`: prepare DEXSeq and run
        [rMATS](http://rnaseq-mats.sourceforge.net/) and
        [MAJIQ](https://majiq.biociphers.org).
      - `drseq_simulation.R`: run DrSeq and DEXSeq, analyze simulation
        results, including all the statistics and plots reported in the
        paper.
      - majiq/: contain two files needed for running MAJIQ.
      - dexseq/: contain two files needed for DEXSeq preparation.

To reproduce the ENCODE results (available at (DOI):
[10.5281/zenodo.3779037](https://doi.org/10.5281/zenodo.3779037)):

1.  Download the processed bam files (shRNA-seq and eCLIP-seq) from
    [ENCODE](https://www.encodeproject.org) portal.
2.  Download transcriptome quantification of
    [TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)
    and [GTEx](https://commonfund.nih.gov/gtex) projects from
    [Xena](https://xena.ucsc.edu).
3.  Run `xena.R`, `encode_surf_one.R` (for each RBP), and
    `encode_surf_summary.R` in order.

To reproduce the simulation results:

1.  Download the processed bam files (Homo sapiens) from ArrayExpress
    dataset
    [E-MTAB-3766](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3766/files/processed/).
2.  Run `other_simulation.sh` and `drseq_simulation.R` in order.

## Contact

Fan Chen (<fan.chen@wisc.edu>) or Sunduz Keles (<keles@stat.wisc.edu>)

## Reference

Chen F and Keles S. “SURF: integrative analysis of a compendium of
RNA-seq and CLIP-seq datasets highlights complex governing of
alternative transcriptional regulation by RNA-binding proteins.”
