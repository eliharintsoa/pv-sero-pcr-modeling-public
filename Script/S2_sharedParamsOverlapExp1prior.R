################################################################################
# SHARED PARAMETERS: global vs village (cluster) model, EXP(1)-PRIOR
#
#
# Output: paper/tables/shared_params_global_vs_cluster_4c_rhoage_exp1prior.csv
#         paper/tables/shared_params_global_vs_cluster_overlap_flags_exp1prior.csv
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})


antigens <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")

# ============================================================================
#shared_params_global_vs_cluster_4c_rhoage_exp1prior.csv
# ============================================================================

ci_fmt <- function(x, digits = 3) {
  fmt <- paste0("%.", digits, "f")
  q   <- quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)
  sprintf(paste0(fmt, " [", fmt, ", ", fmt, "]"), q[2], q[1], q[3])
}

extract_shared_params <- function(ps) {
  if (!is.null(ps$rho) && is.matrix(ps$rho)) {
    rho1 <- ps$rho[, 1]; rho2 <- ps$rho[, 2]; rho3 <- ps$rho[, 3]
  } else if (!is.null(ps[["rho[1]"]])) {
    rho1 <- ps[["rho[1]"]]; rho2 <- ps[["rho[2]"]]; rho3 <- ps[["rho[3]"]]
  } else {
    rho1 <- rho2 <- rho3 <- rep(NA_real_, 1000L)
  }
  list(
    alpha = ps$alpha, rho1 = rho1, rho2 = rho2, rho3 = rho3,
    mu_neg = ps$mu_neg, sigma_neg = ps$sigma_neg, skew_neg = ps$skew_neg,
    delta_sero = ps$delta_sero, sigma_pos = ps$sigma_pos, skew_pos = ps$skew_pos
  )
}

dir_global  <- "skew_normal_fixed_r8_1mu_rho_age"                  
dir_cluster <- "skew_normal_fixed_r8_1mu_rho_age_cluster_exp1prior"  # Exp(1) prior

param_names <- c("alpha", "rho[1]", "rho[2]", "rho[3]",
                 "mu_neg", "sigma_neg", "skew_neg",
                 "delta_sero", "sigma_pos", "skew_pos")

field_map <- list(
  "alpha" = "alpha", "rho[1]" = "rho1", "rho[2]" = "rho2", "rho[3]" = "rho3",
  "mu_neg" = "mu_neg", "sigma_neg" = "sigma_neg", "skew_neg" = "skew_neg",
  "delta_sero" = "delta_sero", "sigma_pos" = "sigma_pos", "skew_pos" = "skew_pos"
)

shared_rows <- list()

for (ag in antigens) {
  
  ps_glob <- tryCatch(readRDS(file.path("results", ag, dir_global,  "posterior_samples.rds")), error = function(e) NULL)
  ps_clus <- tryCatch(readRDS(file.path("results", ag, dir_cluster, "posterior_samples.rds")), error = function(e) NULL)
  if (is.null(ps_glob) && is.null(ps_clus)) next
  
  p_glob <- if (!is.null(ps_glob)) extract_shared_params(ps_glob) else NULL
  p_clus <- if (!is.null(ps_clus)) extract_shared_params(ps_clus) else NULL
  
  for (pname in param_names) {
    fld <- field_map[[pname]]
    v_g <- if (!is.null(p_glob) && !is.null(p_glob[[fld]])) p_glob[[fld]] else NA_real_
    v_c <- if (!is.null(p_clus) && !is.null(p_clus[[fld]])) p_clus[[fld]] else NA_real_
    
    global_str  <- if (all(is.na(v_g))) NA_character_ else ci_fmt(v_g)
    cluster_str <- if (all(is.na(v_c))) NA_character_ else ci_fmt(v_c)
    
    overlap <- NA_real_
    if (!all(is.na(v_g)) && !all(is.na(v_c))) {
      n_samp  <- min(length(v_g), length(v_c))
      overlap <- mean(v_g[1:n_samp] > v_c[1:n_samp])
      overlap <- round(min(overlap, 1 - overlap), 3)
    }
    
    shared_rows[[length(shared_rows) + 1]] <- data.frame(
      antigen = ag, parameter = pname,
      global_model = global_str, cluster_model = cluster_str,
      `P(global>cluster)` = overlap, check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

shared_df <- bind_rows(shared_rows)
dir.create("paper/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(shared_df, "paper/tables/shared_params_global_vs_cluster_4c_rhoage_exp1prior.csv", row.names = FALSE)

