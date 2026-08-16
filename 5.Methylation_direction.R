
######## DIRECTION OF MQTL METHYLATION EFFECT ########

# DO THE MQTL HITS INCREASE OR DECREASE METHYLATION??

library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggrepel)


# LOAD AD GWAS SUBSET
# Subset of SNPs that were mqtl hits to avoid trying to load whole Bellenguez
bellenguez_subset<-read.csv(" ") 

# LOAD LONGEVITY GWAS  SUBSET
# Subset of SNPs that were mqtl hits to avoid trying to load whole GWAS
longevity_subset<-read.csv(" ") 

# LOAD ENTIRE HANNON BRAIN MQTL DATASET (CLEANED)
hannon_mqtls <-read_excel(" ")




########################################
# Add direction column to mqtl df
# Hannon et al., 2016 reported strong concordance in direction of methylation between fetal and adult brain
# Fetal brain regression coefficients therefore used as a representative measure
direction <- function(x)
{case_when(
  x > 0 ~ "increase",
  x < 0 ~ "decrease",
  TRUE  ~ "no effect")}
hannon_mqtls <- hannon_mqtls %>%
  mutate(Fetal_direction = direction(Regression_coefficient_F),
         PFC_direction   = direction(Regression_coefficient_PFC),
         CER_direction   = direction(Regression_coefficient_CER),
         STR_direction   = direction(Regression_coefficient_STR))
table(hannon_mqtls$Fetal_direction)
########################################




##### ALZHEIMER'S DISEASE MQTLS #####

# join AD subset of hits and mqtl df
mqtlhits_merged <- inner_join(bellenguez_subset,hannon_mqtls,by = c("variant_id" = "SNP ID"))
mqtlhits_merged <- mqtlhits_merged %>%
  mutate(flip = effect_allele != Allele_F, beta_harmonised = ifelse(flip, -beta, beta)) #harmonise effect allele



# PLOT 
# Direction of DNA methylation effect vs association with AD

mqtlhits_plot <- ggplot(mqtlhits_merged,
                        aes(x     = Regression_coefficient_F,
                            y     = -log10(p_value), 
                            colour = Fetal_direction,
                            label  = variant_id)) +
  geom_hline(yintercept = -log10(5e-08), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(5e-07), linetype = "dashed", colour = "grey80") +
  geom_hline(yintercept = -log10(5e-06), linetype = "dashed", colour = "grey80") +
  geom_hline(yintercept = -log10(5e-05), linetype = "dashed", colour = "grey80") +
  geom_point(size = 2.5, alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = mqtlhits_merged %>% filter(variant_id %in% key_snps),
    aes(label = variant_id),
    size   = 2.5,
    max.overlaps = 50,
    min.segment.length = Inf,
    box.padding        = 0.2,   
    point.padding      = 0.5, 
    force              = 0.8     
  ) +
  scale_colour_manual(values = c(
    "increase" = "#E74C3C",
    "decrease" = "#3498DB")) +
  labs(
    title   = "AD mQTL hits",
    x       = "mQTL regression coefficient",
    y       = "GWAS association (- log10 p-value)",
    colour  = "DNA methylation"
  ) +
  scale_y_continuous(
    breaks = seq(0, 30, by = 5),
    limits = c(0, 30), expand = c(0, 0)) +
  scale_x_continuous(
    breaks = seq(-0.4, 0.4, by =0.1),
    limits = c(-0.4, 0.4)) +
  theme_bw()
print(mqtlhits_plot)

ggsave("mqtlhits_plot.png", 
       plot = mqtlhits_plot,
       width = 6,
       height = 4,  
       units = "in", dpi = 300)

# Top AD mQTL hits
key_snps <- c("rs5167", "rs16979595", "rs8111069", "rs204468", "rs6743470", "rs199533", "rs7097656", "rs3763312", "rs9394766", "rs1881194", "rs199501")





# PLOT
# Where in the gene are AD mqtl hits? plot distribution of AD brain mQTL hits in gene by CpG position

unique(hannon_mqtls$`Relation to gene`)
simplify_region <- function(x) {case_when(
  is.na(x) ~ NA_character_,
  grepl("TSS200", x) ~ "TSS200",
  grepl("TSS1500", x) ~ "TSS1500",
  grepl("1stExon", x) ~ "1stExon",
  grepl("5'UTR", x) ~ "5'UTR",
  grepl("Body", x) ~ "Body",
  grepl("3'UTR", x) ~ "3'UTR",
  TRUE ~ NA_character_)}

mqtlhits_merged <- mqtlhits_merged %>% mutate(relation_simple = sapply(`Relation to gene`, simplify_region))
table(mqtlhits_merged$relation_simple)
mqtlhits_merged_plot <- mqtlhits_merged %>%
  filter(!is.na(relation_simple),
         Fetal_direction != "no effect") %>%
  count(relation_simple, Fetal_direction) %>%
  tidyr::complete(
    relation_simple,
    Fetal_direction = c("increase", "decrease"),
    fill = list(n = 0))


all_region_plot <- mqtlhits_merged_plot %>%
  mutate(relation_simple = factor(
    relation_simple,
    levels = c("TSS200", "TSS1500", "5'UTR",
               "1stExon", "Body", "3'UTR", "IGR"))) %>%
  ggplot(aes(
    x = relation_simple,
    y = n,
    fill = Fetal_direction)) +
  geom_col(
    position = "dodge", alpha = 0.85) +
  scale_fill_manual(
    values = c(
      "increase" = "#E74C3C",
      "decrease" = "#3498DB")) +
  labs(
    title = "CpG location of all AD-associated mQTL hits",
    x = "CpG relation to gene",
    y = "Number of mQTL hits",
    fill = "DNA methylation") +
  scale_y_continuous(
    breaks = seq(0, 50, by = 5),
    limits = c(0, 50),
    expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic()
print(all_region_plot)

ggsave("all_region_plot.png", 
       plot = all_region_plot,
       width = 6,
       height = 4,  
       units = "in", dpi = 300)

mqtlhits_merged_plot %>%
  tidyr::pivot_wider(
    names_from = Fetal_direction,
    values_from = n,
    values_fill = 0)


########################################

##### Longevity MQTLS #####

# join Longevity subset of hits and mqtl df
long_mqtlhits_merged <- inner_join(longevity_subset,hannon_mqtls,by = c("variant_id" = "SNP ID"))
long_mqtlhits_merged <- long_mqtlhits_merged %>%
  mutate(flip = effect_allele != Allele_F, beta_harmonised = ifelse(flip, -beta, beta)) #harmonise effect allele



# PLOT 
# Direction of DNA methylation effect vs association with Longevity
long_mqtlhits_plot <- ggplot(long_mqtlhits_merged,
                             aes(x     = Regression_coefficient_F,
                                 y     = -log10(p_value), 
                                 colour = Fetal_direction,
                                 label  = variant_id)) +
  ggrepel::geom_text_repel(
    data               = long_mqtlhits_merged %>% filter(variant_id %in% key_long_snps),
    aes(label          = variant_id),
    size               = 2.5,
    min.segment.length = Inf,
    max.overlaps = 50
  ) +
  geom_hline(yintercept = -log10(5e-08), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(5e-07), linetype = "dashed", colour = "grey80") +
  geom_hline(yintercept = -log10(5e-06), linetype = "dashed", colour = "grey80") +
  geom_hline(yintercept = -log10(5e-05), linetype = "dashed", colour = "grey80") +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = c(
    "increase" = "#E74C3C",
    "decrease" = "#3498DB"
  )) +
  labs(
    title   = "Logevity mQTL hits",
    x       = "mQTL regression coefficient",
    y       = "GWAS association (- log10 p-value)",
    colour  = "DNA methylation") +
  scale_y_continuous(
    breaks = seq(0, 30, by = 5),
    limits = c(0, 30), expand = c(0, 0)) +
  scale_x_continuous(
    breaks = seq(-0.4, 0.4, by =0.1),
    limits = c(-0.4, 0.4)) +
  theme_bw()
print(long_mqtlhits_plot)

ggsave("long_mqtlhits_plot.png", 
       plot = long_mqtlhits_plot,
       width = 6,
       height = 4,  
       units = "in", dpi = 300)

# Top longevity mQTL hits
key_long_snps <- c("rs3795437", "rs13099", "rs11130218","rs516246", "rs503279", "rs7845800", "rs130076", "rs1049256", "rs1548807")





# PLOT
# Where in the gene are longevity mqtl hits? plot distribution of longevity brain mQTL hits in gene by CpG position

long_mqtlhits_merged <- long_mqtlhits_merged %>%
  mutate(
    relation_simple = sapply(`Relation to gene`, simplify_region))
table(mqtlhits_merged$relation_simple)
long_mqtlhits_merged_plot <- long_mqtlhits_merged %>%
  filter(!is.na(relation_simple),
         Fetal_direction != "no effect") %>%
  count(relation_simple, Fetal_direction) %>%
  tidyr::complete(
    relation_simple,
    Fetal_direction = c("increase", "decrease"),
    fill = list(n = 0))


long_all_region_plot <- long_mqtlhits_merged_plot %>%
  mutate(relation_simple = factor(
    relation_simple,
    levels = c("TSS200", "TSS1500", "5'UTR",
               "1stExon", "Body", "3'UTR", "IGR"))) %>%
  ggplot(aes(
    x = relation_simple,
    y = n,
    fill = Fetal_direction)) +
  geom_col(
    position = "dodge", alpha = 0.85) +
  scale_fill_manual(
    values = c(
      "increase" = "#E74C3C",
      "decrease" = "#3498DB")) +
  labs(
    title = "CpG location of all AD-associated mQTL hits",
    x = "CpG relation to gene",
    y = "Number of mQTL hits",
    fill = "DNA methylation") +
  scale_y_continuous(
    breaks = seq(0, 50, by = 5),
    limits = c(0, 50),
    expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic()
print(long_all_region_plot)

ggsave("long_all_region_plot.png", 
       plot = long_all_region_plot,
       width = 6,
       height = 4,  
       units = "in", dpi = 300)

long_mqtlhits_merged_plot %>%
  tidyr::pivot_wider(
    names_from = Fetal_direction,
    values_from = n,
    values_fill = 0)
