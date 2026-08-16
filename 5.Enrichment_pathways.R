library(data.table)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(ggh4x)

######### Functional Enrichment Plot #########

# To investigate functional enrichment of mQTLs
# Download results from G:Profiler pathway enrichment analysis and plot graph

# LOAD GO PATHWAY RESULTS
ADmqtl_pathway<-read.csv("/AD_output/pathway_results.csv") 

# CALCULATE GENE RATIO (precision)
ADmqtl_pathway$gene_ratio<-ADmqtl_pathway$intersection_size/ADmqtl_pathway$query_size
colnames(ADmqtl_pathway)

# Order brain regions and type of GO analysis
source_order <- c("GO:CC", "GO:MF", "GO:BP", "KEGG", "MIRNA")
region_order <- c("FETAL", "PFC", "CER", "STR")

ADmqtl_pathway <- ADmqtl_pathway %>%
  filter(!is.na(source), !is.na(region)) %>% 
  mutate(region = factor(region, levels = region_order),
         clumped_status = factor(clumped_status, levels = c("CLUMPED", "NONCLUMPED")),
   source = factor(source, levels = intersect(source_order, unique(source)))
  ) %>% filter(!is.na(region), !is.na(source))


# Order pathway terms: grouped by source, then by highest gene ratio first within each source 
# (most significant G ratio  = highest on y axis within each source panel)
term_order <- ADmqtl_pathway %>%
  group_by(source, term_name) %>%
  summarise(max_ratio = max(gene_ratio, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = factor(source,levels = intersect(source_order,unique(ADmqtl_pathway$source)))) %>%
  arrange(source,
    desc(max_ratio)
  ) %>%
  pull(term_name)

# str_wrap makes different terms into dif blocks
# (insert line breaks into long term names so they wrap at 35 characters instead of running off plot)
ADmqtl_pathway <- ADmqtl_pathway %>%
  mutate(term_name = factor(str_wrap(term_name, 35),
          levels = rev(str_wrap(term_order, 35)))) #highest gene ratio first - reverse order - largest to smallest gene ratio

# PLOT
# add filter to keep only cellular components
pathways_plot <- ggplot(
  data = ADmqtl_pathway %>% filter(source %in% c("GO:CC", "KEGG") |
        (source == "GO:MF" & term_name == "zinc ion binding")),
  aes(x = gene_ratio, y = term_name)) +
  geom_point(
    aes(colour = negative_log10_of_adjusted_p_value,
        shape  = clumped_status,
        size = intersection_size), alpha = 0.85) +
  scale_size_continuous(
    name   = "Gene count",
    range  = c(4, 10),       
    breaks = c(5, 10, 15)) +
  ggh4x::facet_grid2(
    rows      = vars(source),
    cols      = vars(region),
    scales    = "free_y",
    space     = "free_y",         
  ) +
  scale_colour_gradient( #gradient for -log10 p value
    low  = "skyblue",
    high = "royalblue",
    name = expression(-log[10]("adj. p-value"))) +
  scale_shape_manual(
    name   = " ",
    values = c("CLUMPED" = 16, "NONCLUMPED" = 15), #16 codes circles, 15 codes squares
    labels = c("Clumped SNPs", "Non-Clumped SNPs"),
    guide = guide_legend (override.aes = list (size = 6, colour = "royalblue"))) +
  
  scale_x_continuous(
    limits = c(0, 0.5),
    expand = expansion(add = c(0, 0.05)) 
  ) +
  
  labs(x = "Gene Ratio", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    # source strip styling
    strip.text.y       = element_text(face = "bold", size = 12, angle = -90, hjust = 0.5),
    strip.switch.pad.grid = unit(0.5, "cm"),
    strip.text.x       = element_text(face = "bold", size = 12),
    strip.background   = element_rect(fill = "grey92", color = "grey70"),
    panel.spacing.y    = unit(0, "lines"), 
    panel.spacing.x    = unit(0.4, "lines"),
    axis.text.y        = element_text(size = 12),
    axis.text.x        = element_text(size = 10),
    panel.grid.major =  element_line (color = "grey92"),
    legend.position    = "right")

print(pathways_plot)

ggsave("pathways_plot5.png", 
       plot = pathways_plot,
       width = 14,
       height = 7,  
       units = "in", dpi = 300)

