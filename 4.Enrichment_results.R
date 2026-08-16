library(data.table)
library(readxl)
library(dplyr)

setwd("")

# Downlaod Results from HPC for mQTL enrichment
# Repeat with longevity mQTLs

########################################################################################################################

#1. P threshold
#2. Hits (TRUE) = No. of mQTL SNPs which are GWAS significant at each threshold 

#3. qtl rate = proportion of all mQTL SNPs which are significant at that threshold
#4. all rate = proportion of GWAS SNPs which are significant at that threshold
# OR = qtl rate / all rate 

#5. emp p value = significance
#               if I randomly picked the same no. of SNPs from the GWAS (matching for allele frequency),
#               how often would I see this many or more GWAS significant SNPs by chance

########################################################################################################################


# 1.
# ENRICHMENT OF FETAL MQTLS IN AD GWAS
# Load enrichment results downloaded from HPC
ADmqtl_enrichment<-read.csv("Fetal_enrichRESULT_ADmqtl.csv") 

# CLEAN DATA
# Drop first row, so rows have names
rownames(ADmqtl_enrichment) <- ADmqtl_enrichment[[1]]
ADmqtl_enrichment <- ADmqtl_enrichment[, -1]
# Flip the dataframe to make rows columns
ADmqtl_enrichment <- as.data.frame(t(ADmqtl_enrichment))


# CALCULATE ODDS RATIO
# OR = qtl rate / all rate 
ADmqtl_enrichment$OR <- ADmqtl_enrichment$qtl.rate/ADmqtl_enrichment$all.rate
# CALCULATE LOG OR
ADmqtl_enrichment$log2_OR <- log2(ADmqtl_enrichment$OR)


# CALCULATE STANDARD ERROR
# SE of the log OR = √1/a + 1/b + 1/c + 1/d based on 4x4 table

#                Sig  Not sig   
# mQTL SNPs       a    b        total mQtls = 2,564     
# non-mQTL SNPs   c    d        total GWAS - total mQtls = 21,070,465 - 2,564 = 21,067,901
    # a = significant mqtl hits
    ADmqtl_enrichment$a <- ADmqtl_enrichment$'true'
    # b = mqtls that are not significant
    ADmqtl_enrichment$b <- 2564 - ADmqtl_enrichment$a 
    # c = significant non-mqtl hits
                  # rate = sig hits / total
                  # so, total sig hits = total GWAS SNPs x all rate
                  # sig non-mQTL hits = (total sig hits) - sig mQTL hits
    ADmqtl_enrichment$c <- (ADmqtl_enrichment$all.rate * 21070465) - ADmqtl_enrichment$a 
    # d = non-mQTLs that are not significant
    ADmqtl_enrichment$d <- 21070465 - ADmqtl_enrichment$a - ADmqtl_enrichment$b - ADmqtl_enrichment$c
ADmqtl_enrichment$SE <- sqrt(
  1/ADmqtl_enrichment$a + 1/ADmqtl_enrichment$b + 1/ADmqtl_enrichment$c + 1/ADmqtl_enrichment$d)


# CALCULATE 95% CI
# remember to convert back from log scale
ADmqtl_enrichment$log2_OR_CI_low <- ADmqtl_enrichment$log2_OR - 1.96 * ADmqtl_enrichment$SE
ADmqtl_enrichment$log2_OR_CI_high <- ADmqtl_enrichment$log2_OR + 1.96 * ADmqtl_enrichment$SE
ADmqtl_enrichment$CI_low <- exp(ADmqtl_enrichment$log2_OR - 1.96 * ADmqtl_enrichment$SE)
ADmqtl_enrichment$CI_high <- exp(ADmqtl_enrichment$log2_OR + 1.96 * ADmqtl_enrichment$SE)

# drop excess columns to make clean table
  ADmqtl_enrichment <- ADmqtl_enrichment %>% select(-a, -b, -c, -d)

  
  
  
  
# 2.
# ENRICHMENT OF PFC MQTLS IN AD GWAS
 # Load enrichment results downloaded from HPC
   PFC_ADmqtl_enrichment<-read.csv("PFC_enrichRESULT_ADmqtl.csv") 
  
 # CLEAN DATA
  rownames(PFC_ADmqtl_enrichment) <- PFC_ADmqtl_enrichment[[1]]
  PFC_ADmqtl_enrichment <- PFC_ADmqtl_enrichment[, -1]
  PFC_ADmqtl_enrichment <- as.data.frame(t(PFC_ADmqtl_enrichment))
  
 # Calculate OR
  PFC_ADmqtl_enrichment$OR <- PFC_ADmqtl_enrichment$qtl.rate/PFC_ADmqtl_enrichment$all.rate
 # Calculate log OR
  PFC_ADmqtl_enrichment$log2_OR <- log2(PFC_ADmqtl_enrichment$OR)

 # CALCULATE STANDARD ERROR
  # SE of the log OR = √1/a + 1/b + 1/c + 1/d based on 4x4 table (see above)
  #    total mQtls = 1,462
  #    total non-mQTL SNPs = (total GWAS - total mQtls) = (21,070,465 - 1,462) = 21,069,003
       # a
       PFC_ADmqtl_enrichment$a <- PFC_ADmqtl_enrichment$'true'
       # b
       PFC_ADmqtl_enrichment$b <- 1462 - PFC_ADmqtl_enrichment$a 
       # c
       PFC_ADmqtl_enrichment$c <- (PFC_ADmqtl_enrichment$all.rate * 21070465) - PFC_ADmqtl_enrichment$a 
       # d
       PFC_ADmqtl_enrichment$d <- 21070465 - PFC_ADmqtl_enrichment$a - PFC_ADmqtl_enrichment$b - PFC_ADmqtl_enrichment$c
  PFC_ADmqtl_enrichment$SE <- sqrt(
    1/PFC_ADmqtl_enrichment$a + 1/PFC_ADmqtl_enrichment$b + 1/PFC_ADmqtl_enrichment$c + 1/PFC_ADmqtl_enrichment$d)
 # CALCULATE 95% CI
   PFC_ADmqtl_enrichment$log2_OR_CI_low <- PFC_ADmqtl_enrichment$log2_OR - 1.96 * PFC_ADmqtl_enrichment$SE
   PFC_ADmqtl_enrichment$log2_OR_CI_high <- PFC_ADmqtl_enrichment$log2_OR + 1.96 * PFC_ADmqtl_enrichment$SE
   PFC_ADmqtl_enrichment$CI_low <- exp(PFC_ADmqtl_enrichment$log2_OR - 1.96 * PFC_ADmqtl_enrichment$SE)
   PFC_ADmqtl_enrichment$CI_high <- exp(PFC_ADmqtl_enrichment$log2_OR + 1.96 * PFC_ADmqtl_enrichment$SE)
 # drop excess columns to make clean table
   PFC_ADmqtl_enrichment <- PFC_ADmqtl_enrichment %>% select(-a, -b, -c, -d)
  
  

 
 
# 3.
# ENRICHMENT OF CER MQTLS IN AD GWAS
  # Load enrichment results downloaded from HPC
    CER_ADmqtl_enrichment<-read.csv("CER_enrichRESULT_ADmqtl.csv") 
  
 # CLEAN DATA
  rownames(CER_ADmqtl_enrichment) <- CER_ADmqtl_enrichment[[1]]
  CER_ADmqtl_enrichment <- CER_ADmqtl_enrichment[, -1]
  CER_ADmqtl_enrichment <- as.data.frame(t(CER_ADmqtl_enrichment))
  
 # Calculate OR
  CER_ADmqtl_enrichment$OR <- CER_ADmqtl_enrichment$qtl.rate/CER_ADmqtl_enrichment$all.rate
  # Calculate log 2 OR
  CER_ADmqtl_enrichment$log2_OR <- log2(CER_ADmqtl_enrichment$OR)
  
 # CALCULATE STANDARD ERROR
  # SE of the log OR = √1/a + 1/b + 1/c + 1/d based on 4x4 table (see above)
  #    total mQtls = 1,245
  #    total non-mQTL SNPs = (total GWAS - total mQtls) = (21,070,465 - 1,246) = 21,069,219
       # a
       CER_ADmqtl_enrichment$a <- CER_ADmqtl_enrichment$'true'
       # b
       CER_ADmqtl_enrichment$b <- 1246 - CER_ADmqtl_enrichment$a 
       # c
       CER_ADmqtl_enrichment$c <- (CER_ADmqtl_enrichment$all.rate * 21070465) - CER_ADmqtl_enrichment$a 
       # d
       CER_ADmqtl_enrichment$d <- 21070465 - CER_ADmqtl_enrichment$a - CER_ADmqtl_enrichment$b - CER_ADmqtl_enrichment$c
  CER_ADmqtl_enrichment$SE <- sqrt(
    1/CER_ADmqtl_enrichment$a + 1/CER_ADmqtl_enrichment$b + 1/CER_ADmqtl_enrichment$c + 1/CER_ADmqtl_enrichment$d)
 # CALCULATE 95% CI
   CER_ADmqtl_enrichment$log2_OR_CI_low <- CER_ADmqtl_enrichment$log2_OR - 1.96 * CER_ADmqtl_enrichment$SE
   CER_ADmqtl_enrichment$log2_OR_CI_high <- CER_ADmqtl_enrichment$log2_OR + 1.96 * CER_ADmqtl_enrichment$SE
   CER_ADmqtl_enrichment$CI_low <- exp(CER_ADmqtl_enrichment$log2_OR - 1.96 * CER_ADmqtl_enrichment$SE)
   CER_ADmqtl_enrichment$CI_high <- exp(CER_ADmqtl_enrichment$log2_OR + 1.96 * CER_ADmqtl_enrichment$SE)
 # drop excess columns to make clean table
   CER_ADmqtl_enrichment <- CER_ADmqtl_enrichment %>% select(-a, -b, -c, -d)
  
  
  
  
  
# 4.
# ENRICHMENT OF STR MQTLS IN AD GWAS
  # Load enrichment results downloaded from HPC
  STR_ADmqtl_enrichment<-read.csv("STR_enrichRESULT_ADmqtl.csv") 
  
  # CLEAN DATA
  rownames(STR_ADmqtl_enrichment) <- STR_ADmqtl_enrichment[[1]]
  STR_ADmqtl_enrichment <- STR_ADmqtl_enrichment[, -1]
  STR_ADmqtl_enrichment <- as.data.frame(t(STR_ADmqtl_enrichment))
  
  # Calculate OR
  STR_ADmqtl_enrichment$OR <- STR_ADmqtl_enrichment$qtl.rate/STR_ADmqtl_enrichment$all.rate
  # Calculate log 2 OR
  STR_ADmqtl_enrichment$log2_OR <- log2(STR_ADmqtl_enrichment$OR)

  # CALCULATE STANDARD ERROR
  # SE of the log OR = √1/a + 1/b + 1/c + 1/d based on 4x4 table (see above)
  #    total mQtls = 1,350
  #    total non-mQTL SNPs = (total GWAS - total mQtls) = (21,070,465 - 1,350) = 
       # a
       STR_ADmqtl_enrichment$a <- STR_ADmqtl_enrichment$'true'
       # b
       STR_ADmqtl_enrichment$b <- 1350 - STR_ADmqtl_enrichment$a 
       # c
       STR_ADmqtl_enrichment$c <- (STR_ADmqtl_enrichment$all.rate * 21070465) - STR_ADmqtl_enrichment$a 
       # d
       STR_ADmqtl_enrichment$d <- 21070465 - STR_ADmqtl_enrichment$a - STR_ADmqtl_enrichment$b - STR_ADmqtl_enrichment$c
  STR_ADmqtl_enrichment$SE <- sqrt(
    1/STR_ADmqtl_enrichment$a + 1/STR_ADmqtl_enrichment$b + 1/STR_ADmqtl_enrichment$c + 1/STR_ADmqtl_enrichment$d)
  # CALCULATE 95% CI
    STR_ADmqtl_enrichment$log2_OR_CI_low <- STR_ADmqtl_enrichment$log2_OR - 1.96 * STR_ADmqtl_enrichment$SE
    STR_ADmqtl_enrichment$log2_OR_CI_high <- STR_ADmqtl_enrichment$log2_OR + 1.96 * STR_ADmqtl_enrichment$SE
    STR_ADmqtl_enrichment$CI_low <- exp(STR_ADmqtl_enrichment$log2_OR - 1.96 * STR_ADmqtl_enrichment$SE)
    STR_ADmqtl_enrichment$CI_high <- exp(STR_ADmqtl_enrichment$log2_OR + 1.96 * STR_ADmqtl_enrichment$SE)
  # drop excess columns to make clean table
    STR_ADmqtl_enrichment <- STR_ADmqtl_enrichment %>% select(-a, -b, -c, -d)
  
  
  
  

########################################

# PLOT RESULTS

# 1.
# PLOT RESULTS FOR FETAL MQTLS
# Plot OR for each P-value threshold
# ADmqtl_enrichment$pThres <- factor(ADmqtl_enrichment$pThres)
# ADmqtl_enrichment$'TRUE'<- as.numeric(ADmqtl_enrichment$'TRUE')
# ADmqtl_enrichment <- ADmqtl_enrichment %>% rename(hits = `TRUE`)
# ADmqtl_enrichment$significance <- ifelse(ADmqtl_enrichment$emp <0.05, "significant", "not_significant")

# PLOT OF OR
# mqtl_plot <- ggplot(
#  data = ADmqtl_enrichment, 
#  aes(x = pThres, y = OR)) +
#  geom_errorbar(aes( ymin = CI_low, ymax = CI_high ), width = 0.2 ) +
#  geom_point(aes(size = hits)) +
#  scale_size_continuous(breaks = c(5, 10, 15, 20, 25),
#                        limits = c(5, 25),
#                       range = c(1, 4),
#                        name = "mQTL Hits"
#                        ) +
#  geom_hline(yintercept = 1, linetype = "dashed") +
#  labs (x = "GWAS P-value Threshold",
#        y = "Odds Ratio"
#  ) +
#  theme_classic () +
#  theme(
#    axis.title = element_text(size = 10),
#    axis.text = element_text(size = 8)) +
#  scale_y_continuous(
#    breaks = seq(0, 50, by = 5),
#    limits = c(0, 35))
# print(mqtl_plot)


# PLOT OF Log 2 OR
# mqtl_plot2 <- ggplot(
#  data = ADmqtl_enrichment, 
#  aes(x = pThres, y = log2_OR)) +
#  geom_errorbar(aes( ymin = log2_OR_CI_low, ymax = log2_OR_CI_high ), width = 0.2 ) +
#  geom_point(aes(colour = significance, size = hits)) +
#  geom_text(aes(label = hits), hjust = 1.8, size = 3) +
#  scale_colour_manual(values = c("significant" = "skyblue", "not_significant" = "royalblue"),
#                      limits = c("significant", "not_significant"), drop = FALSE,
#                      name = "Emperical p-value", labels = c("significant", "not significant")) +
#  scale_size_continuous(breaks = c(5, 10, 15, 20, 25),
#                        limits = c(1, 25),
#                        range = c(1, 4),
#                        name = "mQTL hits"
#  ) +
#  geom_hline(yintercept = 1, linetype = "dashed") +
#  labs (title = "Enrichment of Fetal mQTLs in AD GWAS",
#        x = "GWAS p-value threshold",
#        y = "Odds Ratio (log2 scale)") +
#  theme_classic () +
#  theme(
#    axis.title = element_text(size = 10),
#    axis.text = element_text(size = 8)) +
#  scale_y_continuous(
#    breaks = seq(0, 50, by = 0.5),
#    limits = c(0.5, 5))
#print(mqtl_plot2)


########################################

# PLOT RESULTS
    
library(ggplot2)

ADmqtl_enrichment$pThres <- factor(ADmqtl_enrichment$pThres)
ADmqtl_enrichment$'true'<- as.numeric(ADmqtl_enrichment$'true')
ADmqtl_enrichment <- ADmqtl_enrichment %>% rename(hits = `true`)
ADmqtl_enrichment$significance <- ifelse(ADmqtl_enrichment$emp <0.05, "significant", "not_significant")

PFC_ADmqtl_enrichment$pThres <- factor(PFC_ADmqtl_enrichment$pThres)
PFC_ADmqtl_enrichment$'true'<- as.numeric(PFC_ADmqtl_enrichment$'true')
PFC_ADmqtl_enrichment <- PFC_ADmqtl_enrichment %>% rename(hits = `true`)
PFC_ADmqtl_enrichment$significance <- ifelse(PFC_ADmqtl_enrichment$emp <0.05, "significant", "not_significant")

CER_ADmqtl_enrichment$pThres <- factor(CER_ADmqtl_enrichment$pThres)
CER_ADmqtl_enrichment$'true'<- as.numeric(CER_ADmqtl_enrichment$'true')
CER_ADmqtl_enrichment <- CER_ADmqtl_enrichment %>% rename(hits = `true`)
CER_ADmqtl_enrichment$significance <- ifelse(CER_ADmqtl_enrichment$emp <0.05, "significant", "not_significant")

STR_ADmqtl_enrichment$pThres <- factor(STR_ADmqtl_enrichment$pThres)
STR_ADmqtl_enrichment$'TRUE'<- as.numeric(STR_ADmqtl_enrichment$'true')
STR_ADmqtl_enrichment <- STR_ADmqtl_enrichment %>% rename(hits = `true`)
STR_ADmqtl_enrichment$significance <- ifelse(STR_ADmqtl_enrichment$emp <0.05, "significant", "not_significant")

# CREATE COMBINED DF
# facet_wrap() for comparing plots

ADmqtl_enrichment$region <- "Fetal"
PFC_ADmqtl_enrichment$region <- "PFC"
CER_ADmqtl_enrichment$region <- "CER"
STR_ADmqtl_enrichment$region <- "STR"

combined_df <- dplyr::bind_rows(
  ADmqtl_enrichment,
  PFC_ADmqtl_enrichment,
  CER_ADmqtl_enrichment,
  STR_ADmqtl_enrichment)

combined_df$region <- factor(combined_df$region, levels = c("Fetal", "PFC", "CER", "STR"))
combined_df$pThres <- factor(combined_df$pThres, levels = c("5e-05", "5e-06", "5e-07", "5e-08"))

combined_plot <- ggplot(
  data = combined_df, 
  aes(x = pThres, y = log2_OR)) +
  geom_errorbar(aes(ymin = log2_OR_CI_low, ymax = log2_OR_CI_high ), width = 0.2 ) +
  geom_point(aes(colour = significance, size = hits)) +
  geom_text(aes(label = hits), hjust = 1.8, size = 3) +
  scale_colour_manual(values = c("significant" = "skyblue", "not_significant" = "royalblue"),
                      limits = c("significant", "not_significant"), drop = FALSE,
                      name = "Emperical p-value", labels = c("significant", "not significant")) +
  scale_size_continuous(breaks = c(5, 10, 15, 20, 25),
                        limits = c(1, 25),
                        range = c(1, 4),
                        name = "mQTL hits"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  labs (x = "GWAS p-value threshold",
        y = "Odds Ratio (log2 scale)"
  ) +
  theme_classic () +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"), strip.placement = "outside",strip.background = element_blank())+
  scale_y_continuous(
    breaks = seq(0, 50, by = 0.5),
    limits = c(0.5, 5)) +
    facet_wrap(~region,  nrow = 1, strip.position ="bottom")
print(combined_plot)

ggsave("combined_plot.png", 
       plot = combined_plot,
       width = 10,
       height = 4,  
       units = "in", dpi = 300)
