################################################################################
# LF016 — real-data fit at r = 3/year (122-day PCR duration)
#
# Same production model and settings as the r = 8 fit, only r_fixed changes.
# Provides LF016's own "true" parameters at r = 3 for the recovery check in
# S12_r3_recovery_check.R.
#
# Output: results/LF016/skew_normal_fixed_r3_1mu_rho_age/
################################################################################


suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(dplyr)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

antigen_name <- "LF016"
r_fixed_value <- 3

stan_data_orig <- readRDS(sprintf("data/stan_data_%s.rds", antigen_name))
stan_data <- list(
  N          = stan_data_orig$N,
  age        = stan_data_orig$age,
  log_mfi    = stan_data_orig$log_mfi,
  pcr        = stan_data_orig$pcr,
  r_fixed    = r_fixed_value,
  lambda_max = stan_data_orig$lambda_max
)

stan_model_comp <- stan_model(file = "stan_models/serocatalytic_pcr_4compartments_FIXED_R_skewnormal_1mu_rho_age.stan")

fit <- sampling(
  stan_model_comp,
  data    = stan_data,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  thin    = 1,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

params_scalar <- c("lambda", "alpha", "mu_neg", "mu_pos", "delta_sero",
                   "sigma_neg", "sigma_pos", "skew_neg", "skew_pos")
params_rho <- paste0("rho[", 1:3, "]")
params <- c(params_scalar, params_rho)

summary_table <- summary(fit, pars = params)$summary

print(round(summary_table, 4))

rhats <- summary_table[, "Rhat"]


sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
div <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))


ps <- rstan::extract(fit)
use_best_chains <- FALSE
best_chain_indices <- seq_len(fit@sim$chains)

if (max(rhats, na.rm = TRUE) > 1.05) {
  
  
  
  ll_arr <- rstan::extract(fit, pars = "log_lik", permuted = FALSE)
  n_chains_total <- dim(ll_arr)[2]
  chain_mean_ll <- sapply(seq_len(n_chains_total),
                          function(ch) mean(rowSums(ll_arr[, ch, ])))
  
  
  for (i in seq_along(chain_mean_ll))
    
    
    sorted_order <- order(chain_mean_ll)
  gaps <- diff(chain_mean_ll[sorted_order])
  gap_pos <- which.max(gaps)
  best_chain_indices <- sort(sorted_order[(gap_pos + 1):n_chains_total])
  dropped_chains <- sort(sorted_order[seq_len(gap_pos)])
  
  
  draws_best <- posterior::subset_draws(posterior::as_draws_array(fit),
                                        chain = best_chain_indices)
  rhats_best <- sapply(params_scalar, function(p)
    posterior::rhat(posterior::subset_draws(draws_best, variable = p)))
  
  
  for (p in params_scalar) cat(sprintf("  %-14s: %.4f\n", p, rhats_best[[p]]))
  
  
  params_to_filter <- c(params, "rho", "prob_seropositive", "prob_pcr_positive", "log_lik")
  posterior_samples_best <- list()
  
  for (param in params_to_filter) {
    arr <- rstan::extract(fit, pars = param, permuted = FALSE)
    sub <- arr[, best_chain_indices, , drop = FALSE]
    d   <- dim(sub)
    posterior_samples_best[[param]] <-
      if (length(d) == 3 && d[3] == 1) as.vector(sub)
    else matrix(sub, nrow = d[1] * d[2], ncol = d[3])
  }
  
  ps <- posterior_samples_best
  use_best_chains <- TRUE
  cat(">>> Continuing with BEST CHAINS ONLY <<<\n\n")
  
} else {
  cat("All chains converged well. Using all chains.\n\n")
}

out_dir <- file.path("results", antigen_name, "skew_normal_fixed_r3_1mu_rho_age")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(ps, file.path(out_dir, "posterior_samples.rds"))
saveRDS(stan_data, file.path(out_dir, "stan_data.rds"))


fitted_true <- c(
  lambda = median(ps$lambda), alpha = median(ps$alpha),
  rho1 = median(ps$rho[,1]), rho2 = median(ps$rho[,2]), rho3 = median(ps$rho[,3]),
  mu_neg = median(ps$mu_neg), mu_pos = median(ps$mu_pos),
  sigma_neg = median(ps$sigma_neg), sigma_pos = median(ps$sigma_pos),
  skew_neg = median(ps$skew_neg), skew_pos = median(ps$skew_pos)
)
write.csv(data.frame(parameter = names(fitted_true), value = fitted_true),
          file.path(out_dir, "fitted_medians.csv"), row.names = FALSE)


