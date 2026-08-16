install.packages("coloc")
library(coloc)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggrepel)

# ATTEMPT AT RUNNING COLOCALISATION ANALYSIS

# type = cc (case control) for GWAS, quant (quantitative) for mqtls

#H0: neither trait has a genetic association in the region
#H1: only trait 1 (GWAS) has a genetic association in the region
#H2: only trait 2 (mqtl) has a genetic association in the region
#H3: both traits are associated, but with different causal SNPs
#H4: both traits are associated and share a single causal SNP

# LOAD BELLENGUEZ SUBSET
# Subset of SNPs that were mqtl hits to avoid trying to load whole Bellenguez
bellenguez_subset<-read.csv("/Bellenguez_hits_HannonFetal_nonclumped.csv") 
colnames(bellenguez_subset)

# LOAD ENTIRE HANNON BRAIN MQTL DATASET (CLEANED)
hannon_mqtls <-read_excel(" ")
colnames(hannon_mqtls)

################################################################################

######## SET UP GWAS ########
gwas_coloc <- list(
beta   = bellenguez_subset$beta,
varbeta = bellenguez_subset$standard_error^2, #variance = standard error squared
N      =  788989,  # ( CASES PLUS CONTROLS - 111,326 clinically diagnosed/‘proxy’ AD cases and 677,663 controls )
type   = "cc",          
snp    = bellenguez_subset$variant_id)

######## SET UP mQTLs ########
# Calculate varbeta (standard error) by working backwards
hannon_mqtls <- hannon_mqtls %>%mutate(
    z_score_F = qnorm(P_value_F / 2, lower.tail = FALSE),
    SE_F      = abs(Regression_coefficient_F / z_score_F),
    varbeta_F = SE_F^2)
#QUANTITATIVE ANALYSIS REQUIRES MAF
#gwas dataset has MAF
merged <- inner_join(bellenguez_subset,hannon_mqtls, by = c("variant_id" = "SNP ID"))
merged$MAF <- pmin(merged$effect_allele_frequency, 1 - merged$effect_allele_frequency)
# SET UP MQTL DATASET
mqtl_coloc <- list(
  beta   = merged$Regression_coefficient_F, #beta is named regression coefficient in Hannon df
  varbeta = merged$varbeta_F,
  N      =  173,  # ( 173 fetal brain samples (94 male, 79 female) used for DNA methylation and SNP profiling )
  MAF = merged$MAF,
  type   = "quant",          
  snp    = merged$"variant_id")

################################################################################

######## RESULTS ########

# Test using entire genome
result <- coloc.abf(gwas_coloc, mqtl_coloc)
result$summary
# nsnps    PP.H0.abf    PP.H1.abf    PP.H2.abf    PP.H3.abf    PP.H4.abf 
# 1.960000e+02 3.291595e-29 1.249232e-07 2.634895e-22 9.999998e-01 1.011954e-07 
# "H3 = both datasets have signal but driven by different casual SNP"

#SCREEN FOR COLOCALISED SITES
coloc_site=result$results %>% filter(SNP.PP.H4 > 0.8)
# rs5167
# SNP.PP.H4 = 0.9348706
# MAPS TO APOC4 (Chr19)

# WANTED TO PERFORM COLOC IN APOC4 REGION
# standard coloc window is 500kb
# CANNOT ACTUALLY PERFORM COLOC ANALYSIS AND SUSIE AS DO NOT HAVE NON SIG MQTLS IN HANNON DATASET
# Instead had a further look at other SNPs in this region

################################################################################







######## INSPECTING APOE/APOC4 LOCUS ########

# rs5167 = chr19:44945208



#LOAD IN ENTIRE AREA AROUND APOE
chr_19<-read.csv("/Bellenguez_chr19.csv") 
apoc_locus <-chr_19 %>%    
  filter(chromosome == 19,
         base_pair_location >= 44945208 - 250000	,
         base_pair_location <= 44945208 + 250000	)
cat("GWAS SNPs in APOC region:", nrow(apoc_locus), "\n")  # 8155 variants

apoc4_merged <- left_join(apoc_locus,hannon_mqtls,by = c("variant_id" = "SNP ID"))
apoc4_merged <- apoc4_merged %>% mutate(
  methylation_status = case_when(
      Fetal_direction == "increase" ~ "increase",
      Fetal_direction == "decrease" ~ "decrease",
      TRUE ~ "no mqtl association"))

# TOP 5 AD GENES IN THIS REGION = "rs747519137", "rs117310449", "rs144261139", 	"rs139644294", "rs537741299"




# SUBSET AREA AROUND APOE FOR MQTL HITS ONLY
apoc_hits <- bellenguez_subset %>%    
  filter(chromosome == 19,
         base_pair_location >= 44945208 - 250000	,
         base_pair_location <= 44945208 + 250000	)
cat("GWAS SNPs in APOC region:", nrow(apoc_hits), "\n")

# 5 SNPs
# "rs747519137", "rs117310449", "rs144261139", 	"rs139644294", "rs537741299"


key_apoc_snps <- c("rs5167", "rs16979595", "rs7253458", "rs8111069", "rs204468", "rs747519137", "rs117310449", "rs144261139", 	"rs139644294", "rs537741299")

# Plot regional position against GWAS association
regional_apoc_plot <- ggplot(apoc4_merged,
                        aes(x      = base_pair_location,
                            y      = -log10(p_value),
                            colour = methylation_status)) +
  geom_point(
    data  = apoc4_merged %>% filter(methylation_status == "no mqtl association"),
    size  = 1.5, alpha = 0.4) +
  geom_point(
    data  = apoc4_merged %>% filter(methylation_status != "no mqtl association"),
    size  = 2.5, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = apoc4_merged %>% filter(variant_id %in% key_apoc_snps),
    aes(label = variant_id),
    colour = "black",
    size   = 2.5) +
  geom_hline(yintercept = -log10(5e-08), linetype = "dashed", colour = "black") +
  scale_colour_manual(values = c(
    "increase" = "#E74C3C",   # red
    "decrease" = "#3498DB",   # blue
    "no mqtl association"   = "grey80")) +
  labs(
    title  = "APOE locus",
    x      = "Chr19 position",
    y      = "-log10 p-value",
    colour = "DNA methylation") +
  scale_y_continuous(
    breaks = seq(0, 120, by = 20),
    limits = c(0, 120), expand = c(0, 0)) +
  theme_bw()
print(regional_apoc_plot)

ggsave("regional_apoc_plot.png", 
       plot = regional_apoc_plot,
       width = 10,
       height = 6,  
       units = "in", dpi = 300)




################################################################################
######## INSPECTING MAPT LOCUS ########


#LOAD IN ENTIRE AREA AROUND MAPT
chr_17<-read.csv("/Bellenguez_chr17.csv") 
mapt_locus <- chr_17 %>%    
  filter(chromosome == 17,
         base_pair_location >= 45894278  - 500000,
         base_pair_location <= 46028334  + 500000	)
cat("GWAS SNPs in MAPT region:", nrow(mapt_locus), "\n") #6191 variants

mapt_merged <- left_join(mapt_locus, hannon_mqtls, by = c("variant_id" = "SNP ID"))
mapt_merged <- mapt_merged %>% mutate(
    methylation_status = case_when(
      Fetal_direction == "increase" ~ "increase",
      Fetal_direction == "decrease" ~ "decrease",
      TRUE                          ~ "no mqtl association"))




# SUBSET AREA AROUND MAPT FOR MQTL HITS ONLY
mapt_hits <- bellenguez_subset %>%    
filter(chromosome == 17,
       base_pair_location >= 45894278  - 500000,
       base_pair_location <= 46028334  + 500000	)
# 31 SNPs
key_mapt_snps <- c("rs1881194", "rs7350928", "rs17574604", "rs1052587", "rs17574361", "rs17577094", "rs1052553", "rs12185233", "rs1052551", "rs17691610")
key_mapt_snps2 <- c("rs1881194", "rs7350928", "rs17574604", "rs1052587", "rs17574361", "rs17577094", "rs1052553", "rs12185233", "rs1052551", "rs17691610", 
                   "rs17651549", "rs17660464", "rs12373142", "rs17652121", "rs12185235", "rs12185268", "rs12373123", "rs17763596", "rs16940665", "rs1396862", "rs17334797", "rs17426064")



regional_mapt_plot <- ggplot(mapt_merged,
                        aes(x      = base_pair_location,
                            y      = -log10(p_value),
                            colour = methylation_status)) +
  geom_point(
    data  = mapt_merged %>% filter(methylation_status == "no mqtl association"),
    size  = 1.5, alpha = 0.4
  ) +
  geom_point(
    data  = mapt_merged %>% filter(methylation_status != "no mqtl association"),
    size  = 2.5, alpha = 0.9
  ) +
  geom_hline(yintercept = -log10(5e-08), linetype = "dashed", colour = "black") +
  ggrepel::geom_text_repel(
    data = mapt_merged %>% filter(variant_id %in% key_mapt_snps),
    aes(label = variant_id),
    colour = "black",
    size   = 2.5, max.overlaps = 50,
  ) +
  scale_colour_manual(values = c(
    "increase" = "#E74C3C",   # red
    "decrease" = "#3498DB",   # blue
    "no mqtl association"   = "grey80"
  )) +
  labs(
    title  = "MAPT locus",
    x      = "Chr17 position",
    y      = "-log10 p-value",
    colour = "DNA methylation"
  ) +
  scale_y_continuous(
    breaks = seq(0, 120, by = 20),
    limits = c(0, 120), expand = c(0, 0)) +
  theme_bw()
print(regional_mapt_plot)

ggsave("regional_mapt_plot.png", 
       plot = regional_mapt_plot,
       width = 10,
       height = 6,  
       units = "in", dpi = 300)