# HPC code for enrichment analysis, all-cause GWAS 
library(readxl)
library(foreach)
library(doParallel)
library(parallel)
library(data.table)

# For use in HPC 
# Load in AD GWAS (background)
gwas_all <- fread("/Documents/enrichment_analysis/AD_input/filtered_gwas.csv")
# 21,070,465

# Load in mQTL FOREGROUND 
gwas_qtl <- read.csv("/Documents/enrichment_analysis/AD_input/fetal_mqtl_foreground.csv")
#2,564 out of 2,604 fetal mQTLs are found in the GWAS
#1,462 out of 1,482 PFC mQTLs are found in the GWAS
#1,246 out of 1,266 CER mQTLs are found in the GWAS
#1,350 out of 1,374 STR mQTLs are found in the GWAS
#147,118 out of 148,826 blood mQRLs are found in the GWAS

pThres <- c(5e-5, 5e-6, 5e-7, 5e-8)


### ---- DEFINE MAF BINS CORRECTLY ----

bin_width <- 0.02
breaks <- seq(0, 1, by = bin_width)

# Assign bins ONCE using cut()
gwas_all$bin  <- cut(gwas_all$MAF,
                     breaks = breaks,
                     include.lowest = TRUE,
                     right = TRUE)

gwas_qtl$bin  <- cut(gwas_qtl$MAF,
                     breaks = breaks,
                     include.lowest = TRUE,
                     right = TRUE)


# Count how many mQTL SNPs per bin
freq_bins <- table(gwas_qtl$bin)

# Remove empty bins to avoid indexing confusion
freq_bins <- freq_bins[freq_bins > 0]

# Create GWAS SNP pools per bin
bin_pools <- split(seq_len(nrow(gwas_all)), gwas_all$bin)

# Keep only bins used by QTLs
bin_pools <- bin_pools[names(freq_bins)]

# Sanity check
stopifnot(sum(freq_bins) == nrow(gwas_qtl))


### ---- PERMUTATION FUNCTION  ----

permutation <- function(pThres, freq_bins, bin_pools, gwas_all) {
  
  sim_indices <- integer(sum(freq_bins))
  pos <- 1
  
  for (b in names(freq_bins)) {
    
    n_needed <- as.integer(freq_bins[b])
    pool <- bin_pools[[b]]
    
    sel <- sample(pool, n_needed, replace = TRUE)
    
    sim_indices[pos:(pos + n_needed - 1)] <- sel
    pos <- pos + n_needed
  }
  
  gwas_sub <- gwas_all[sim_indices, ]
  
  sapply(pThres, function(p)
    sum(gwas_sub$P < p)
  )
}



################################################################################
# FOR HPC ONLY
# Create parallel cluster here, after defining objects above
# Detect number of cores allocated by scheduler, and cap this at 4 to try not to exceed memory
ncores <- min(as.numeric(Sys.getenv("NSLOTS")), 4)

if (is.na(ncores) || ncores <1) {
  ncores <- 1 
}

cat("Cores allocated:", ncores, "\n")

# Prevent thread oversubscription inside workers
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

cl <- makeCluster(ncores, type = "PSOCK")
registerDoParallel(cl)

clusterExport(cl, c("permutation","pThres","freq_bins","bin_pools","gwas_all"))
################################################################################


# Set no. of permutations
# In HPC run 10,000 in 100 batches, when running in local R use a smaller no.

# Run permutations 
nPerms <- 1000
batch_size <- 100
n_batches <- nPerms / batch_size

perms.tmp <- matrix(NA_integer_, nrow = length(pThres), ncol = nPerms)
col_idx <- 1

for (b in 1:n_batches) {
  message("Running batch ", b, " of ", n_batches)
  
  batch_res <- foreach(
    i = 1:batch_size,
    .combine  = cbind,
    .packages = c("stats"),
    .errorhandling = "stop", 
    .options.snow = list(preschedule = FALSE)
  ) %dopar% {
    permutation(pThres, freq_bins, bin_pools, gwas_all)
  }
  
  perms.tmp[, col_idx:(col_idx + batch_size - 1)] <- batch_res
  col_idx <- col_idx + batch_size
}


dim(perms.tmp)

length(pThres)

true<-vector(length = 4)
emp<-vector(length = 4)
all<-vector(length = 4)

for(k in 1:length(pThres)){
  
  all[k]<-length(which(gwas_all$P < pThres[k]))
  true[k]<-length(which(gwas_qtl$P < pThres[k]))
  emp[k] <- (sum(perms.tmp[k,] >= true[k], na.rm = TRUE) + 1) /
    (sum(!is.na(perms.tmp[k,])) + 1)
}

qtl.rate<-true/nrow(gwas_qtl)
all.rate<-all/nrow(gwas_all)

print(true)
print(nrow(gwas_qtl))
print(qtl.rate)


# Save Results
save(perms.tmp, true, pThres, file = "permutation_results_ADmqtl.RData")
write.csv(rbind(pThres, true,qtl.rate, all.rate, emp), "Fetal_enrichRESULT_ADmqtl.csv")

stopCluster(cl)


################################################################################

# To identify signficant mQTL hits
# Run nonclumped analysis with HPC.R script and 10 permutations in local R, 
### ---- EXTRACT mQTL SNPs PER GWAS THRESHOLD ----

# Create a list to store SNP data frames per threshold
mqtl_thresholds <- lapply(seq_along(pThres), function(k) {
  
  df <- gwas_qtl[gwas_qtl$P < pThres[k], ]
  
  # Add a column indicating which threshold this came from
  df$threshold <- pThres[k]
  
  return(df)
})
sapply(pThres, function(x)
  sum(gwas_qtl$P < x)
)

# Name the list elements for clarity
names(mqtl_thresholds) <- paste0("p<", pThres)

# Combine into one single data frame (long format)
mqtl_snps_df <- do.call(rbind, mqtl_thresholds)

# Optional: reset row names
rownames(mqtl_snps_df) <- NULL

# Quick check
head(mqtl_snps_df)

# Save list of SNPs
write.csv(mqtl_snps_df, "nonclumped_AD_mQTL_hits.csv")
