#!/usr/bin/env Rscript
################################################################################
# SENSITIVITY ANALYSIS: FIXED r, POPULATION-LEVEL 4C-rho-age MODEL (EBP)
#
# Conducted on a single reference antigen (EBP): the mechanism this tests -
# the closed-form dependence of age-specific PCR prevalence on lambda and r
# alone, independent of alpha and rho - is determined by age and PCR status,
# identical across all eight antigen-specific fits, rather than by
# antigen-specific serological parameters.
# Output: results/EBP/sensitivity_r_rhoage/
#   sensitivity_summary.csv, loo_summary.csv, pcr_curves_by_r.rds
#   elpd_vs_r.pdf, parameters_vs_r.pdf, pcr_prevalence_curves_by_r.pdf
################################################################################

rm(list = ls())
# Run from the repository root (Pv_sero_pcr_modeling/)

antigen_name <- "EBP"
r_values     <- c(0.5, 1, 2, 3, 4, 6, 8)

suppressPackageStartupMessages({
  library(rstan)
  library(loo)
  library(dplyr)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("compartment_model.R")   # calculate_compartments_age_varying()

################################################################################
# LOAD DATA
################################################################################

data_file <- sprintf("Data/stan_data_%s.rds", antigen_name)
stan_data_orig <- readRDS(data_file)

################################################################################
# COMPILE MODEL ONCE (original, unmodified - exponential(10.0) prior on lambda,
# already confirmed prior-insensitive at population level)
################################################################################

stan_file <- "Stan/4comp_rho_age.stan"
stan_model_compiled <- stan_model(file = stan_file)

################################################################################
# FIT FOR EACH r — extract only small summaries, discard fit/posterior/loo objects
################################################################################

summary_rows <- list()
loo_rows     <- list()
age_seq      <- seq(2, 95, by = 0.5)
pcr_curve_list <- list()

for (r_val in r_values) {
  
  stan_data <- list(
    N          = stan_data_orig$N,
    age        = stan_data_orig$age,
    log_mfi    = stan_data_orig$log_mfi,
    pcr        = stan_data_orig$pcr,
    r_fixed    = r_val,
    lambda_max = stan_data_orig$lambda_max
  )
  
  start_time <- Sys.time()
  fit <- sampling(
    stan_model_compiled,
    data    = stan_data,
    chains  = 4,
    iter    = 2000,
    warmup  = 1000,
    thin    = 1,
    cores   = 4,
    seed    = 42,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  elapsed <- as.numeric(Sys.time() - start_time, units = "mins")
  
  post <- rstan::extract(fit, pars = c("lambda", "alpha", "rho"))
  summary_table <- summary(fit, pars = c("lambda", "alpha", "rho"))$summary
  
  log_lik    <- loo::extract_log_lik(fit, parameter_name = "log_lik")
  loo_result <- loo(log_lik, cores = 4)
  
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  rhats <- summary_table[, "Rhat"]
  ess   <- summary_table[, "n_eff"]
  
  rho_med <- apply(post$rho, 2, median)
  
  summary_rows[[as.character(r_val)]] <- data.frame(
    r_fixed = r_val,
    pcr_duration_days = round(365 / r_val),
    lambda_median = median(post$lambda),
    lambda_lo = quantile(post$lambda, 0.025), lambda_hi = quantile(post$lambda, 0.975),
    alpha_median = median(post$alpha),
    alpha_lo = quantile(post$alpha, 0.025), alpha_hi = quantile(post$alpha, 0.975),
    rho1_median = rho_med[1], rho2_median = rho_med[2], rho3_median = rho_med[3],
    hl1_months = log(2)/rho_med[1]*12, hl2_months = log(2)/rho_med[2]*12, hl3_months = log(2)/rho_med[3]*12,
    elapsed_min = elapsed, divergences = divergences,
    max_rhat = max(rhats, na.rm = TRUE), min_ess = min(ess, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  loo_rows[[as.character(r_val)]] <- data.frame(
    r_fixed = r_val,
    elpd_loo = loo_result$estimates["elpd_loo","Estimate"],
    elpd_se  = loo_result$estimates["elpd_loo","SE"],
    p_loo    = loo_result$estimates["p_loo","Estimate"],
    looic    = loo_result$estimates["looic","Estimate"],
    n_pareto_k_high = sum(loo_result$diagnostics$pareto_k > 0.7)
  )
  
  # PCR-prevalence-by-age curve at posterior-median params, for the plot only
  lam <- median(post$lambda); alp <- median(post$alpha)
  pcr_curve_list[[as.character(r_val)]] <- sapply(age_seq, function(a) {
    S <- calculate_compartments_age_varying(a, lam, alp, rho_med, r_val)
    S["N1"] + S["P1"]
  })
  
  # discard
  rm(fit, post, log_lik, loo_result, sampler_params, summary_table)
  gc(verbose = FALSE)
}

################################################################################
# SAVE — small files 
################################################################################

out_dir <- file.path("results", antigen_name, "sensitivity_r_rhoage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- do.call(rbind, summary_rows); rownames(summary_df) <- NULL
loo_df     <- do.call(rbind, loo_rows);     rownames(loo_df) <- NULL
loo_df$elpd_diff_vs_r8 <- loo_df$elpd_loo - loo_df$elpd_loo[loo_df$r_fixed == 8]

print(summary_df, row.names = FALSE, digits = 4)
print(loo_df, row.names = FALSE, digits = 4)

write.csv(summary_df, file.path(out_dir, "sensitivity_summary.csv"), row.names = FALSE)
write.csv(loo_df,     file.path(out_dir, "loo_summary.csv"),         row.names = FALSE)
saveRDS(pcr_curve_list, file.path(out_dir, "pcr_curves_by_r.rds"))  # small: 7 x length(age_seq) vectors

################################################################################
# PLOTS (small PDFs only)
################################################################################

pdf(file.path(out_dir, "elpd_vs_r.pdf"), width = 8, height = 6)
plot(loo_df$r_fixed, loo_df$elpd_loo, type = "b", pch = 16, col = "darkblue", lwd = 2, cex = 1.5,
     xlab = "r_fixed (per year)", ylab = "ELPD (LOO-CV)",
     main = paste("Model fit (LOO-CV) vs r_fixed -", antigen_name))
arrows(loo_df$r_fixed, loo_df$elpd_loo - 2*loo_df$elpd_se,
       loo_df$r_fixed, loo_df$elpd_loo + 2*loo_df$elpd_se, angle = 90, code = 3, length = 0.1, col = "darkblue")
axis(3, at = summary_df$r_fixed, labels = paste0(summary_df$pcr_duration_days, "d"), col.axis = "grey50", cex.axis = 0.8)
mtext("PCR duration (days)", side = 3, line = 2.5, col = "grey50")
grid()
best_idx <- which.max(loo_df$elpd_loo)
points(loo_df$r_fixed[best_idx], loo_df$elpd_loo[best_idx], pch = 1, cex = 3, col = "red", lwd = 2)
legend("bottomright", "Best model", pch = 1, col = "red", pt.cex = 2, pt.lwd = 2)
dev.off()


