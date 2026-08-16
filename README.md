# VBLR9_MSc_Clinical_Neuroscience
Integrating Genomic and Epigenomic Analyses to Investigate Mechanisms Underlying Late Onset Alzheimer’s Disease and Ageing

This repository contains all of the analysis scripts used for my MSc project (2025/2026)

1. ## mQTL clumping
   Clean mQTL datasets
   Clump mQTLs to account for linkage disequilibrium (LD) using PLINK2 (MHC region clumped separately)
   
2. ## AD and longevity preprocessing
   Clean GWAS datasets
   Create mQTL foreground using clumped mQTLs for use in HPC enrichment analyses (subset of all mQTL SNPs also present in GWAS dataset).
   
3. ## AD and longevity enrichment analysis HPC scripts
   Designed to be run on a HPC cluster via shell script (mqtl_enrichment.sh)
   
4. ## AD and longevity enrichment results
   Download and process results from HPC (calculate OR and 95% CI, plot enrichment graphs)
   
5. ## Create plots
   - Functional enrichment of mQTL hits (results downloaded from G:Profiler)
   - Visualisation of mQTL hits using heatmaps
   - Plot of methylation effects against GWAS association and plot of distribution of significant mQTL hits by CpG position 
   - Colocalisation screening and follow up investigation of APOC4 and MAPT loci

## Shell script: used to run R enrichment analyses on HPC



Note: working directories and file paths within each script will need to be edited to match local environment before running.
