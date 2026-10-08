################################################################################
# BUILD STAN_DATA FILES FROM Data/data_cleaned.csv
#
# Generates the population-level and cluster-stratified Stan input files that
# every model-fitting script (02-06) and the r-sensitivity script (11) read
# via readRDS(). Run this once before fitting any model.
#
# lambda_max / r_max are fixed prior-support bounds (not derived from data).
#
# Output: Data/stan_data_<antigen>.rds              (population-level)
#         Data/stan_data_<antigen>_stratified.rds    (cluster-level)
################################################################################

# Run from the repository root (Pv_sero_pcr_modeling/)

antigens <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")

build_stan_data <- function(antigen_name, stratified = FALSE,
                            data_file = "Data/data_cleaned.csv") {
  d <- read.csv(data_file)
  
  mfi_col <- paste0("log_mfi_", antigen_name)
  
  out <- list(
    N       = nrow(d),
    age     = d$agen,
    log_mfi = d[[mfi_col]],
    pcr     = as.integer(d$pcr_positive)
  )
  
  if (stratified) {
    commune_f <- factor(d$commune)
    out$N_commune     <- length(levels(commune_f))
    out$commune       <- as.integer(commune_f)
    out$commune_names <- levels(commune_f)
    out$N_cluster      <- length(unique(d$cluster))
    out$cluster        <- d$cluster
    out$cluster_names  <- as.character(sort(unique(d$cluster)))
  }
  
  out$lambda_max <- 6
  out$r_max      <- 100
  
  out
}

for (ag in antigens) {
  saveRDS(build_stan_data(ag, stratified = FALSE),
          file.path("Data", sprintf("stan_data_%s.rds", ag)))
  saveRDS(build_stan_data(ag, stratified = TRUE),
          file.path("Data", sprintf("stan_data_%s_stratified.rds", ag)))
}
