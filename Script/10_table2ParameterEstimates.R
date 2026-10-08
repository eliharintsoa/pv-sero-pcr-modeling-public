# =============================================================================
# TABLE 2 — Population-level parameter estimates
#
# Uses posterior_samples_best_chains.rds when available
#
# Outputs:
#   paper/table2_4comp_rho_age_parameters.csv   — Table 2 (4C rho-age model), wide format
#   paper/table2_2comp_parameters.csv           — 2-compartment model, supplementary
#   paper/table2_combined.csv                   — key params, both constant-rho models
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

# Run from the repository root (Pv_sero_pcr_modeling/)
dir.create("paper", showWarnings = FALSE)

ANTIGENS <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")

# Load best chains if available, otherwise fall back to all chains
load_ps <- function(dir) {
  f_best <- file.path(dir, "posterior_samples_best_chains.rds")
  f_all  <- file.path(dir, "posterior_samples.rds")
  if (file.exists(f_best)) {
    ps <- readRDS(f_best); attr(ps, "source") <- "best_chains"; return(ps)
  }
  if (file.exists(f_all)) {
    ps <- readRDS(f_all);  attr(ps, "source") <- "all_chains";  return(ps)
  }
  NULL
}

fmt_ci <- function(draws, digits = 3) {
  fmt_ci_vals(median(draws), quantile(draws, 0.025), quantile(draws, 0.975), digits)
}
fmt_ci_vals <- function(med, lo, hi, digits = 3) {
  sprintf("%.*f (%.*f, %.*f)", digits, med, digits, lo, digits, hi)
}

# =============================================================================
# 4-compartment model
# =============================================================================
params_4c <- c("lambda", "alpha", "rho[1]", "rho[2]", "rho[3]",
               "mu_neg", "mu_pos", "sigma_neg", "sigma_pos",
               "skew_neg", "skew_pos")

tbl_4c <- do.call(rbind, lapply(ANTIGENS, function(ag) {
  ps <- load_ps(file.path("results", ag, "skew_normal_fixed_r8_1mu_rho_age"))
  if (is.null(ps)) { warning("Missing posterior for ", ag); return(NULL) }
  out <- data.frame(Antigen = ag, stringsAsFactors = FALSE)
  for (p in params_4c) {
    # Handle rho[1], rho[2], rho[3]
    if (grepl("^rho\\[([123])\\]$", p)) {
      idx <- as.integer(sub("^rho\\[([123])\\]$", "\\1", p))
      rho_mat <- ps[["rho"]]
      if (is.null(rho_mat)) { out[[p]] <- NA_character_; next }
      draws <- if (is.matrix(rho_mat)) rho_mat[, idx] else rho_mat[idx]
    } else {
      draws <- ps[[p]]
      if (is.null(draws)) { out[[p]] <- NA_character_; next }
      if (is.matrix(draws)) draws <- as.vector(draws)
    }
    out[[p]] <- fmt_ci(draws)
  }
  out
}))

write.csv(tbl_4c, "paper/table2_4comp_rho_age_parameters.csv", row.names = FALSE)

# =============================================================================
# 2-compartment model
# =============================================================================
params_2c_raw <- c("lambda", "rho", "mu_0", "mu_1", "epsilon_0", "epsilon_1",
                   "alpha_0", "alpha_1")

tbl_2c <- do.call(rbind, lapply(ANTIGENS, function(ag) {
  ps <- load_ps(file.path("results", ag, "2compartments_skewnormal"))
  if (is.null(ps)) { warning("Missing posterior for ", ag); return(NULL) }
  out <- data.frame(Antigen = ag, stringsAsFactors = FALSE)
  for (p in params_2c_raw) {
    draws <- ps[[p]]
    if (is.null(draws)) { out[[p]] <- NA_character_; next }
    if (is.matrix(draws)) draws <- as.vector(draws)
    out[[p]] <- fmt_ci(draws)
  }
  # Derived: mu_pos = mu_0 + mu_1
  out[["mu_pos"]] <- fmt_ci(ps$mu_0 + ps$mu_1)
  out
}))

write.csv(tbl_2c, "paper/table2_2comp_parameters.csv", row.names = FALSE)

# =============================================================================
# Combined: key parameters only
# =============================================================================
tbl_key_4c <- do.call(rbind, lapply(ANTIGENS, function(ag) {
  ps <- load_ps(file.path("results", ag, "skew_normal_fixed_r8"))
  if (is.null(ps)) return(NULL)
  data.frame(Model = "4-compartment", Antigen = ag,
             lambda = fmt_ci(ps$lambda),
             alpha  = fmt_ci(ps$alpha),
             rho    = fmt_ci(ps$rho),
             stringsAsFactors = FALSE)
}))

tbl_key_2c <- do.call(rbind, lapply(ANTIGENS, function(ag) {
  ps <- load_ps(file.path("results", ag, "2compartments_skewnormal"))
  if (is.null(ps)) return(NULL)
  data.frame(Model = "2-compartment", Antigen = ag,
             lambda = fmt_ci(ps$lambda),
             alpha  = NA_character_,
             rho    = fmt_ci(ps$rho),
             stringsAsFactors = FALSE)
}))

tbl_combined <- rbind(tbl_key_4c, tbl_key_2c)
tbl_combined <- tbl_combined[order(tbl_combined$Antigen, tbl_combined$Model), ]
write.csv(tbl_combined, "paper/table2_combined.csv", row.names = FALSE)

