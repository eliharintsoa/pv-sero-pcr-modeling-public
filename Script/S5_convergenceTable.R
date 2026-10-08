# ============================================================================
# MCMC CONVERGENCE DIAGNOSTICS — all four models, all eight antigens
#
# R-hat: Gelman-Rubin statistic (target < 1.01). ESS_min: minimum effective
# sample size across all parameters (target > 400). Statistics are computed
# on the well-mixed chains only for the constant-rho models
# (fit_best_chains.rds, same chains used for all parameter estimates and
# credible intervals; 
# Output: paper/tables/convergence_diagnostics.csv
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(rstan)
  library(posterior)
})


dir.create("paper/tables", showWarnings = FALSE, recursive = TRUE)

antigens <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")

model_dirs <- c(
  "2comp_const"   = "2compartments_skewnormal",
  "2comp_rho_age" = "skew_normal_2comp_rho_age",
  "4comp_const"   = "skew_normal_fixed_r8",
  "4comp_rho_age" = "skew_normal_fixed_r8_1mu_rho_age"
)
model_labels <- c(
  "2comp_const"   = "2C-const",
  "2comp_rho_age" = "2C-rho-age",
  "4comp_const"   = "4C-const",
  "4comp_rho_age" = "4C-rho-age (preferred)"
)
const_rho_models <- c("2comp_const", "4comp_const")

fit_path_fn <- list(
  "2comp_const"   = function(ag) sprintf("results/%s/2compartments_skewnormal/fit.rds", ag),
  "2comp_rho_age" = function(ag) sprintf("results/%s/skew_normal_2comp_rho_age/fit.rds", ag),
  "4comp_const"   = function(ag) sprintf("results/%s/skew_normal_fixed_r8/fit_fixed_r.rds", ag),
  "4comp_rho_age" = function(ag) sprintf("results/%s/skew_normal_fixed_r8_1mu_rho_age/fit.rds", ag)
)

# For constant-rho models: load pre-filtered draws_array (fit_best_chains.rds)
# and compute diagnostics via posterior::summarise_draws() — same chains used
# for all parameter estimates and credible intervals.
best_chain_stats <- function(ag, mkey) {
  dir       <- sprintf("results/%s/%s", ag, model_dirs[mkey])
  path_best <- file.path(dir, "fit_best_chains.rds")
  path_fit  <- fit_path_fn[[mkey]](ag)
  
  if (file.exists(path_best)) {
    # Filtered draws_array: use it directly
    draws <- tryCatch(readRDS(path_best), error = function(e) NULL)
    if (is.null(draws)) return(NULL)
  } else if (file.exists(path_fit)) {
    # All chains mixed well — use full stanfit converted to draws_array
    fit <- tryCatch(readRDS(path_fit), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    draws <- tryCatch(posterior::as_draws_array(fit), error = function(e) NULL)
    if (is.null(draws)) return(NULL)
  } else {
    return(NULL)
  }
  # Keep only scalar model parameters (no [i] arrays); exclude lp__
  all_vars   <- posterior::variables(draws)
  par_vars   <- all_vars[!grepl("[", all_vars, fixed = TRUE) & all_vars != "lp__"]
  draws_pars <- posterior::subset_draws(draws, variable = par_vars)
  summ <- tryCatch(
    posterior::summarise_draws(draws_pars,
                               rhat     = posterior::rhat,
                               ess_bulk = posterior::ess_bulk),
    error = function(e) NULL)
  if (is.null(summ)) return(NULL)
  list(
    Rhat_max = round(max(summ$rhat,     na.rm = TRUE), 4),
    ESS_min  = round(min(summ$ess_bulk, na.rm = TRUE))
  )
}

conv_rows <- list()
for (ag in antigens) {
  for (mkey in names(fit_path_fn)) {
    
    if (mkey %in% const_rho_models) {
      # Use pre-filtered best-chain draws — same chains as parameter estimates
      st <- best_chain_stats(ag, mkey)
      if (is.null(st)) {
        conv_rows[[length(conv_rows)+1]] <- data.frame(
          Antigen=ag, Model=model_labels[mkey],
          Rhat_max=NA, ESS_min=NA, stringsAsFactors=FALSE)
      } else {
        conv_rows[[length(conv_rows)+1]] <- data.frame(
          Antigen   = ag, Model = model_labels[mkey],
          Rhat_max  = st$Rhat_max, ESS_min = st$ESS_min,
          stringsAsFactors = FALSE)
      }
    } else {
      path <- fit_path_fn[[mkey]](ag)
      if (!file.exists(path)) {
        conv_rows[[length(conv_rows)+1]] <- data.frame(
          Antigen=ag, Model=model_labels[mkey],
          Rhat_max=NA, ESS_min=NA, stringsAsFactors=FALSE)
        next
      }
      fit <- tryCatch(readRDS(path), error=function(e) NULL)
      if (is.null(fit)) next
      s <- tryCatch(summary(fit)$summary, error=function(e) NULL)
      if (is.null(s)) next
      conv_rows[[length(conv_rows)+1]] <- data.frame(
        Antigen    = ag,
        Model      = model_labels[mkey],
        Rhat_max   = round(max(s[,"Rhat"],  na.rm=TRUE), 4),
        ESS_min    = round(min(s[,"n_eff"], na.rm=TRUE)),
        stringsAsFactors=FALSE
      )
    }
  }
}

conv_df <- bind_rows(conv_rows)
colnames(conv_df) <- c("Antigen","Model","Rhat_max","ESS_min")

print(conv_df, row.names = FALSE)

write.csv(conv_df, "paper/tables/convergence_diagnostics.csv", row.names = FALSE)
