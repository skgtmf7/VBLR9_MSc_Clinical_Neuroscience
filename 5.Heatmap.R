library(ggplot2)
library(dplyr)
library(tidyr)

######## AD MQTL HEATMAP ########

setwd("/AD_output ")

# Load in extracted SNP hits
# fetalhits<-read.csv("AD_mQTL_hits.csv")
# PFChits<-read.csv("PFC_AD_mqtl_hits.csv")
# CERhits<-read.csv("CER_AD_mqtl_hits.csv")
# STRhits<-read.csv("STR_AD_mqtl_hits.csv")

# Load in extracted SNP hits NONCLUMPED
fetalhits<-read.csv("nonclumped_AD_mqtl_hits.csv")
PFChits<-read.csv("PFC_nonclumped_AD_mqtl_hits.csv")
CERhits<-read.csv("CER_nonclumped_AD_mqtl_hits.csv")
STRhits<-read.csv("STR_nonclumped_AD_mqtl_hits.csv")

# Clean dataframes
  rownames(fetalhits) <- fetalhits[[1]]
  fetalhits <- fetalhits[, -1]
  rownames(PFChits) <- PFChits[[1]]
  PFChits <- PFChits[, -1]
  rownames(CERhits) <- CERhits[[1]]
  CERhits <- CERhits[, -1]
  rownames(STRhits) <- STRhits[[1]]
  STRhits <- STRhits[, -1]

  
  
# Load in annotated gene lists
# fetalhits_genes<-read.csv("fetalhits_genelist.csv", skip = 1)
# PFChits_genes<-read.csv("PFChits_genelist.csv", skip = 1)
# CERhits_genes<-read.csv("CERhits_genelist.csv", skip = 1)
# STRhits_genes<-read.csv("STRhits_genelist.csv", skip = 1)

# Load in annotated gene lists NONCLUMPED
fetalhits_genes<-read.csv("nonclumped_fetalhits_genelist.csv", skip = 1)
PFChits_genes<-read.csv("nonclumped_PFC_genelist.csv", skip = 1)
CERhits_genes<-read.csv("nonclumped_CER_genelist.csv", skip = 1)
STRhits_genes<-read.csv("nonclumped_STR_genelist.csv", skip = 1)



# Add gene list to extracted SNP dataframe (matched by rsid)
fetalhits$gene_names <- fetalhits_genes$gene_names[match(fetalhits$SNP, fetalhits_genes$id)]
PFChits$gene_names <- PFChits_genes$gene_names[match(PFChits$SNP, PFChits_genes$id)]
CERhits$gene_names <- CERhits_genes$gene_names[match(CERhits$SNP, CERhits_genes$id)]
STRhits$gene_names <- STRhits_genes$gene_names[match(STRhits$SNP, STRhits_genes$id)]

# Add region column
fetalhits$region <- "Fetal"
PFChits$region   <- "PFC"
CERhits$region   <- "CER"
STRhits$region   <- "STR"
# COMBINE DATAFRAMES
all_hits <- bind_rows(fetalhits, PFChits, CERhits, STRhits)


########################################

# Layout heatmap table
heatmap_df <- all_hits %>%
  mutate(gene_names = gsub("ENSG\\d+", "", gene_names, perl = TRUE)) %>%
  mutate(gene_names = gsub("^[\\s.,;]+|[\\s.,;]+$", "", gene_names, perl = TRUE))
  filter(!is.na(gene_names), gene_names != "", ) %>%
  group_by(gene_names, region) %>%
  summarise(
    threshold = min(as.numeric(as.character(threshold))),
    .groups = "drop")
heatmap_df$threshold <- factor(heatmap_df$threshold,
                               levels = c(5e-05, 5e-06, 5e-07, 5e-08))
heatmap_df$region <- factor(heatmap_df$region,
                            levels = c("Fetal", "PFC", "CER", "STR"))

# PLOT HEATMAP
ADmqtl_heatmap <- ggplot(
  data = heatmap_df,
  aes(x = region, y = gene_names, fill = threshold)) +
  geom_tile(color = "white", width = 1.2, height = 0.9) + #creates hit map boxes, white border
  scale_fill_manual(
    values = c(
      "5e-05" = "#C6E2FF",
      "5e-06" = "#B9D3EE",
      "5e-07" = "#9FB6CD",
      "5e-08" = "#6C7B8B"),
       na.value = "") +
  theme_classic() +
  labs(x = "Region",
       y = "Gene",
       fill = "GWAS p-value threshold") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12)
        )

print(ADmqtl_heatmap)

ggsave("ADmqtl_heatmap.png", 
       plot = ADmqtl_heatmap,
       width = 10,
       height = 20,  
       units = "in", dpi = 300)


################################################################################

######## LONGEVITY MQTL HEATMAP ########

setwd("/Longevity_output")

# Load in extracted SNP hits NONCLUMPED
long_fetalhits<-read.csv("nonclumped_longevity_mqtl_hits.csv")
long_PFChits<-read.csv("PFC_nonclumped_longevity_mqtl_hits.csv")
long_CERhits<-read.csv("CER_nonclumped_longevity_mqtl_hits.csv")
long_STRhits<-read.csv("STR_nonclumped_longevity_mqtl_hits.csv")
# Clean dataframes
rownames(long_fetalhits) <- long_fetalhits[[1]]
long_fetalhits <- long_fetalhits[, -1]
rownames(long_PFChits) <- long_PFChits[[1]]
long_PFChits <- long_PFChits[, -1]
rownames(long_CERhits) <- long_CERhits[[1]]
long_CERhits <- long_CERhits[, -1]
rownames(long_STRhits) <- long_STRhits[[1]]
long_STRhits <- long_STRhits[, -1]


# Load in annotated gene lists NONCLUMPED
long_fetalhits_genes<-read.csv("nonclumped_longevity_fetal_genelist.csv", skip = 1)
long_PFChits_genes<-read.csv("nonclumped_longevity_PFC_genelist.csv", skip = 1)
long_CERhits_genes<-read.csv("nonclumped_longevity_CER_genelist.csv", skip = 1)
long_STRhits_genes<-read.csv("nonclumped_longevity_STR_genelist.csv", skip = 1)


# Add gene list to extracted SNP dataframe (matched by rs id)
long_fetalhits$gene_names <- long_fetalhits_genes$gene_names[match(long_fetalhits$SNP, long_fetalhits_genes$id)]
long_PFChits$gene_names <- long_PFChits_genes$gene_names[match(long_PFChits$SNP, long_PFChits_genes$id)]
long_CERhits$gene_names <- long_CERhits_genes$gene_names[match(long_CERhits$SNP, long_CERhits_genes$id)]
long_STRhits$gene_names <- long_STRhits_genes$gene_names[match(long_STRhits$SNP, long_STRhits_genes$id)]

# Add region column
long_fetalhits$region <- "Fetal"
long_PFChits$region   <- "PFC"
long_CERhits$region   <- "CER"
long_STRhits$region   <- "STR"
long_BLOODhits$region   <- "Blood"
# COMBINE DATAFRAMES
all_LONGEVITY_hits <- bind_rows(long_fetalhits, long_PFChits, long_CERhits,long_STRhits)


########################################

# Layout heatmap table
longevity_heatmap_df <- all_LONGEVITY_hits %>% 
  mutate(gene_names = gsub("ENSG\\d+", "", gene_names, perl = TRUE)) %>%
  mutate(gene_names = gsub("^[\\s.,;]+|[\\s.,;]+$", "", gene_names, perl = TRUE))
  filter(!is.na(gene_names), gene_names != "") %>%
  group_by(gene_names, region) %>%
  summarise(
    threshold = min(threshold),
    .groups = "drop")
longevity_heatmap_df$threshold <- factor(longevity_heatmap_df$threshold,
                               levels = c(5e-05, 5e-06, 5e-07, 5e-08))
longevity_heatmap_df$region <- factor(longevity_heatmap_df$region,
                            levels = c("Fetal", "PFC", "CER", "STR"))

# PLOT HEATMAP
longevity_mqtl_heatmap <- ggplot(
  data = longevity_heatmap_df,
  aes(x = region, y = gene_names, fill = threshold)) +
  geom_tile(color = "white", height = 0.9) +
  scale_fill_manual(
    values = c(
      "5e-05" = "#FFC1C4",
      "5e-06" = "#FFC1D8",
      "5e-07" = "#CD919E",
      "5e-08" = "#8B6969"),
    drop = FALSE, na.value = "") +
  theme_classic() +
  labs(x = "Region",
       y = "Gene",
       fill = "GWAS p-value threshold") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12)
  )

print(longevity_mqtl_heatmap)

ggsave("longevity_mqtl_heatmap.png", 
       plot = longevity_mqtl_heatmap,
       width = 8,
       height = 12,  
       units = "in", dpi = 300)


