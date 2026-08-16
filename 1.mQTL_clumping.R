######### clump mQTL datsets into semi-independent SNPs using PLINK2 #########

library(readxl)
library(data.table)
library(bigsnpr)
library(bigstatsr)
library(rlang)
library(data.table)
setwd(" ")


# Each mQTL dataset clumped seperately (fetal brain, PFC, STR, CER, blood)
# Load in brain mQTL dataset (Hannon et al., 2016)
non_clumped_mQTL_dataset <- read_excel("Hannon_NonImputed_mQTLs.xlsx")
# Load in blood mQTL dataset (Villicana et al., 2023)
non_clumped_mQTL_dataset <- read.csv("blood_mqtl_clump_input.csv")

# Load PLINK
setwd("/Documents/MSc_code/GWAS/bin/")
plink2 <- "/Documents/MSc_code/GWAS/bin/plink2"
system(paste(plink2, "--version"))

# Load in reference genome
# see details of this, is the 1000 Genomes project Florian Prive, non-related common SNPs
# see github for code by which this was generated, 'common' means in HapMap3 and UKBB 
# could use other reference panel for clumping if appropriate eg. subset for specific popualtion
# or select another dataset if appropriate 
download_1000G("ref_panel")

# Prefix of 1000G reference panel (NO extensions)
bfile <- "ref_panel/1000G_phase3_common_norel"

# Output prefix
out_prefix <- "clumping_results"

########################################
######## CLEAN DATA ########

# Set SNP ID column 'SNP'
setnames(non_clumped_mQTL_dataset,
         "rsid",
         "SNP")

# Set P value column 'P'
# For Brain mQTL datasets, there are P value columns for each brain region:
# (P_value_F, P_value_PFC, P_value_STR, P_value_CER)
# Start with fetal brain region
# Edit with adult regions
setnames(non_clumped_mQTL_dataset,
         "P_value_F",
         "P")
# Edit for blood mQTLs
#setnames(non_clumped_mQTL_dataset,
#         "p.value",
#         "P")
colnames(non_clumped_mQTL_dataset)


#### Remove duplicated SNPs in mQTL df and keep only most significant p value row per duplicate 
# Look if there are duplicated SNPs
dup_SNPs <- unique(non_clumped_mQTL_dataset[["SNP"]][duplicated(non_clumped_mQTL_dataset[["SNP"]])])
length(dup_SNPs) #check how many duplicates
# Order based on p value (smallest to largest)
non_clumped_clean<-non_clumped_mQTL_dataset[order(non_clumped_mQTL_dataset$SNP, non_clumped_mQTL_dataset$P), ]
# CHECK (by finding an id with duplicates and the row with smallest p value should be first)
dup_test <- dup_SNPs[1]
non_clumped_clean[non_clumped_clean$SNP == dup_test, c("SNP", "P")]
# Remove duplicates
non_clumped_clean<-non_clumped_clean[!duplicated(non_clumped_clean$SNP), ]
# CHECK
sum(duplicated(non_clumped_clean$SNP)) #Should be 0

#### Remove missing p values
non_clumped_clean <- non_clumped_clean[!non_clumped_clean$P %in% c("nt", "ns"),]
clump_input <- data.frame(SNP=non_clumped_clean$SNP, P=non_clumped_clean$P)
clump_input <- clump_input[!is.na(clump_input$SNP) & !is.na(clump_input$P), ]
head(clump_input)

fwrite(clump_input, "clump_mQTLs_input.txt", sep = "\t")



########################################

# Genome wide clumping (MHC region clumped separately)

# Specify MHC region 
mhc_range <- data.frame(
  CHR = 6,
  START = as.integer(25000000),
  END = as.integer(35000000),
  LABEL = "MHC")
write.table(mhc_range,
            file = "mhc_region.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = " ")
readLines("mhc_region.txt")


# clump genome except MHC region
cmd_genome <- paste(
  plink2,
  "--bfile ref_panel/1000G_phase3_common_norel",
  "--clump clump_mQTLs_input.txt",
  "--clump-p1 1",
  "--clump-p2 1",
  "--clump-r2 0.25",
  "--clump-kb 250",
  "--chr 1-22",                         # keep all autosomes
  "--exclude range mhc_region.txt",     # remove ONLY MHC region
  "--out", paste0(out_prefix, "_fetal_genome"))
system(cmd_genome)


# MHC extract region for clumping 
cmd_mhc_extract <- paste(
  plink2,
  "--bfile ref_panel/1000G_phase3_common_norel",
  "--chr 6",
  "--from-bp 25000000",
  "--to-bp 35000000",
  "--make-bed",
  "--out", paste0(out_prefix, "_mhc_region"))
system(cmd_mhc_extract)


# clump MHC regions 
cmd_mhc_clump <- paste(
  plink2,
  "--bfile", paste0(out_prefix, "_mhc_region"),
  "--clump clump_mQTLs_input.txt",
  "--clump-p1 1",
  "--clump-p2 1",
  "--clump-r2 0.25",
  "--clump-kb 10000",
  "--out", paste0(out_prefix, "_fetal_mhc"))
system(cmd_mhc_clump)


cmd_mhc_clump<- paste(
  plink2,
  "--bfile", paste0(out_prefix, "_mhc_region"),
  "--clump clump_mQTLs_input.txt",
  "--clump-p1 1",
  "--clump-p2 1",
  "--clump-r2 0.25",
  "--clump-kb 10000",
  "--out", paste0(out_prefix, "1"))
system(cmd_mhc_clump)

########################################

# Combine the MHC clumped region with the rest of the data

get_index_snps <- function(file) {
  dat <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  dat$SNP}

clumped <- fread("clumping_results_fetal_genome.clumps") 
clumped_mhc <- fread("clumping_results_fetal_mhc.clumps")

# double check that these are made correctly ie no overlapping and SNP numbers add up as expected
overlap_snps <- intersect(clumped, clumped_mhc)

snps_genome <- clumped$SNP
snps_mhc <- clumped_mhc$SNP
overlap_snps <- intersect(snps_genome, snps_mhc)
length(overlap_snps) # Yes 0 overlapping SNPs

final_snps <- unique(c(snps_genome, snps_mhc))

length(final_snps) 
# Fetal brain = 2,604 mQTLs
# Adult PFC = 1,482 mQTLs
# Adult CER = 1,266 mQTLs
# Adult STR = 1,374 mQTLs
# Blood = 148,826 mQTLs

fwrite(data.table(SNP = final_snps),
       "fetlal_clumped_mQTL_list.txt",
       col.names = FALSE)

# Make subset of initial mQTL data frame containing only these SNPs to get other columns back
# (careful to subset data frame which does not include duplicates)
clumped_mqtl <- non_clumped_clean[non_clumped_clean$SNP %in% final_snps,]
dim(clumped_mqtl)

fwrite(
  clumped_mqtl,
  "/Documents/MSc_code/enrichment_analysis/fetal_clumped_mQTLs.txt",
  col.names = TRUE,
  sep ="\t")

# Ready for downstream analysis :) 


