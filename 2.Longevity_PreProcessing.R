######### Pre processing steps #########
# Clean Longevity GWAS dataframe
# Create mQTL foreground (subset of GWAS SNPs that are also present in the mQTL dataset)

library(data.table)
library(readxl)
setwd(" ")

# Load Longevity GWAS dataset (Pilling et al., 2017)
# GWAS data frame has been substted to include only relevant columns (SNP, freq, p value), to use less computational power in future analyses
gwas<-read.table("Pilling.tsv", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
# 11,521,815 SNPs
colnames(gwas)

# Load mQTL datasets
# These have been clumped to create quasi-indpendent SNPs
# (start with fetal, repeat for each mqtl dataset)
mQTLs<-fread("fetal_clumped_mQTLs.txt") 
# Fetal brain = 2,604 mQTLs
# Adult PFC = 1,482 mQTLs
# Adult CER = 1,266 mQTLs
# Adult STR = 1,374 mQTLs
# Blood = 148,826 mQTLs
colnames(mQTLs)

################################################################################

# Clean Longevity GWAS dataframe

gwas_2<-as.data.table(gwas) #took too long in base R, so converted to data table

#1. Remove SNPs with unknown rsid, p-value or MAF
gwas_2<-gwas_2[!is.na(gwas_2$variant_id) &
                 !is.na(gwas_2$p_value) &
                 !is.na(gwas_2$effect_allele_frequency)]

#2. Remove duplicates based on best p value
sum(duplicated(gwas_2$variant_id))
# order dataset based on p value (smallest to largest)
gwas_2<-gwas_2[order(variant_id, p_value)]
# CHECK by finding variant id with duplicates and the row with smallest p value should be first
head(gwas_2)
gwas_2[gwas_2$variant_id == gwas_2$variant_id[duplicated(gwas_2$variant_id)][2], ]
# remove duplicates
gwas_3<-gwas_2[!duplicated(gwas_2$variant_id), ]
# CHECK duplicates have been removed
sum(duplicated(gwas_3$variant_id)) #Should be 0

#3. Make sure the MAF is annotated
# standardize across GWAS by keeping whichever is the smallest [effect allele] OR [1 - effect allele] to get MAF
gwas_3$MAF <- pmin(gwas_3$effect_allele_frequency, 1 - gwas_3$effect_allele_frequency)

#4. Create new subset with rsid, p value and MAF to make cleaner without the effect allele column
colnames(gwas_3)
gwas_4 <- data.frame(SNP = gwas_3$variant_id,
                     P = gwas_3$p_value,
                     MAF = gwas_3$MAF)

################################################################################

# CREATE MQTL FOREGROUND
# subset of the GWAS data frame which contains only SNPs present in the mQTL dataset
overlap_SNPs<-(intersect(mQTLs$SNP, gwas_4$SNP))
mqtl_foreground<-gwas_4[gwas_4$SNP %in% overlap_SNPs, ]
#2,472 out of 2,604 fetal mQTLs are found in the GWAS
#1,415 out of 1,482 PFC mQTLs are found in the GWAS
#1,209 out of 1,266 CER mQTLs are found in the GWAS
#1,297 out of 1,374 STR mQTLs are found in the GWAS
#141,460 out of 148,826 blood mQTLs are found in the GWAS

# save filtered gwas dataframe
write.csv(gwas_4, "longevity_filtered_gwas.csv", row.names = FALSE)

# save filtered mqtl foreground
write.csv(mqtl_foreground, "fetal_longevity_mqtl_foreground.csv", row.names = FALSE)

# Ready for mQTL enrichment analysis on HPC :)

################################################################################

# NON CLUMPED MQTL FOREGROUNDS
# To create heatmaps and investigate g:Profiler pathway analysis etc.

# Load (clean) non clumped mQTLs
non_clumped_mQTLs <- read.table("clump_mQTLs_input.txt")
colnames(non_clumped_mQTLs)

# Make mQTL foreground of non clumped mQTLs
overlap_SNPs<-(intersect(non_clumped_mQTLs$SNP, gwas_4$SNP))
nonclumped_mqtl_foreground<-gwas_4[gwas_4$SNP %in% overlap_SNPs, ]
#9,828 out of 10,305 fetal mQTLs are found in the GWAS
#5,203 out of 5,455 PFC mQTLs are found in the GWAS
#4,573 out of 4,767 CER mQTLs are found in the GWAS
#4,755 out of 5,004 STR mQTLs are found in the GWAS

#save filtered (non clumped) mqtl foreground
write.csv(nonclumped_mqtl_foreground, "fetal_nonclumped_longevity_mqtl_foreground.csv", row.names = FALSE)
