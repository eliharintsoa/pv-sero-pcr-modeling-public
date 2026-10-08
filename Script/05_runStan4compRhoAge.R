################################################################################
# RUN STAN — 4-COMPARTMENT MODEL
#   Fixed r | Skew-normal MFI | Age-varying rho | 1 mu_pos
#
# Model: stan_models/serocatalytic_pcr_4compartments_FIXED_R_skewnormal_1mu_rho_age.stan
#
# Age groups: 0–5 (rho[1]) | 5–18 (rho[2]) | >18 (rho[3])
#
#
# Output: results/{antigen}/skew_normal_fixed_r{r}_1mu_rho_age/
################################################################################

rm(list = ls())
# Run from the repository root (Pv_sero_pcr_modeling/)

if (!exists("antigen_name"))  antigen_name  <- "RBP2b.P87"
if (!exists("r_fixed_value")) r_fixed_value <- 8

################################################################################
# 1. LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(rstan)
  library(bayesplot)
  library(loo)
  library(posterior)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Chain log-likelihood diagnostic plot
plot_chain_loglik <- function(fit, antigen, model_tag) {
  ll_arr   <- rstan::extract(fit, pars = "log_lik", permuted = FALSE)
  n_iter   <- dim(ll_arr)[1]
  n_chains <- dim(ll_arr)[2]
  palette  <- c("#2166ac", "#d6604d", "#f4a582", "#4393c3",
                "#92c5de", "#b2182b", "#fddbc7", "#053061")
  colours  <- palette[seq_len(n_chains)]
  labels   <- paste0("Chain ", seq_len(n_chains))
  
  df_ll <- data.frame(
    iteration = rep(seq_len(n_iter), n_chains),
    chain     = factor(rep(labels, each = n_iter), levels = labels),
    total_ll  = as.vector(apply(ll_arr, c(1, 2), sum))
  )
  
  ggplot(df_ll, aes(x = iteration, y = total_ll, colour = chain)) +
    geom_line(alpha = 0.8, linewidth = 0.5) +
    scale_colour_manual(values = setNames(colours, labels)) +
    labs(
      title  = sprintf("Chain log-likelihood — %s | %s", antigen, model_tag),
      x      = "Iteration",
      y      = "Total log-likelihood",
      colour = "Chain"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title       = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

################################################################################
# 2. LOAD DATA
################################################################################

data_file <- sprintf("Data/stan_data_%s.rds", antigen_name)

stan_data_orig <- readRDS(data_file)

stan_data <- list(
  N          = stan_data_orig$N,
  age        = stan_data_orig$age,
  log_mfi    = stan_data_orig$log_mfi,
  pcr        = stan_data_orig$pcr,
  r_fixed    = r_fixed_value,
  lambda_max = if (!is.null(stan_data_orig$lambda_max)) stan_data_orig$lambda_max else 2.0
)

################################################################################
# 3. COMPILE AND RUN
################################################################################

stan_file <- "Stan/4comp_rho_age.stan"

start_compile   <- Sys.time()
stan_model_comp <- stan_model(file = stan_file)

log_mfi_min   <- min(stan_data$log_mfi)
log_mfi_max   <- max(stan_data$log_mfi)
log_mfi_range <- log_mfi_max - log_mfi_min

init_fn <- function() {
  list(
    lambda     = 0.3,
    alpha      = 0.5,
    rho        = c(0.5, 0.2, 0.1),
    mu_neg     = 2.0,
    delta_sero = 4.0,
    sigma_neg  = 1.0,
    sigma_pos  = 1.5,
    skew_neg   = 1.5,
    skew_pos   = -5.0
  )
}

start_sampling <- Sys.time()

fit <- sampling(
  stan_model_comp,
  data    = stan_data,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  thin    = 1,
  cores   = 4,
  seed    = 42,
  # init    = init_fn,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

sampling_time <- as.numeric(Sys.time() - start_sampling, units = "mins")

################################################################################
# 4. CONVERGENCE DIAGNOSTICS
################################################################################

params_scalar <- c("lambda", "alpha",
                   "mu_neg", "mu_pos", "delta_sero",
                   "sigma_neg", "sigma_pos", "skew_neg", "skew_pos")
params_rho    <- paste0("rho[", 1:3, "]")
params        <- c(params_scalar, params_rho)

summary_table <- summary(fit, pars = params)$summary
print(round(summary_table, 4))

rhats    <- summary_table[, "Rhat"]
ess_bulk <- summary_table[, "n_eff"]

if (max(rhats, na.rm = TRUE) >= 1.05) {
  print(names(rhats)[rhats > 1.05])
}

sampler_params     <- get_sampler_params(fit, inc_warmup = FALSE)
divergences        <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
treedepth_hits     <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 15)))

################################################################################
# 4b. CHAIN SELECTION (if convergence is poor)
################################################################################

posterior_samples  <- rstan::extract(fit)
use_best_chains    <- FALSE
best_chain_indices <- seq_len(fit@sim$chains)

if (max(rhats, na.rm = TRUE) > 1.05) {
  
  ll_arr         <- rstan::extract(fit, pars = "log_lik", permuted = FALSE)
  n_chains_total <- dim(ll_arr)[2]
  chain_mean_ll  <- sapply(seq_len(n_chains_total),
                           function(ch) mean(rowSums(ll_arr[, ch, ])))
  
  sorted_order       <- order(chain_mean_ll)
  gaps               <- diff(chain_mean_ll[sorted_order])
  gap_pos            <- which.max(gaps)
  best_chain_indices <- sort(sorted_order[(gap_pos + 1):n_chains_total])
  dropped_chains     <- sort(sorted_order[seq_len(gap_pos)])
  ll_diff            <- mean(chain_mean_ll[best_chain_indices]) -
    mean(chain_mean_ll[dropped_chains])
  
  within_range <- if (length(best_chain_indices) > 1)
    diff(range(chain_mean_ll[best_chain_indices])) else 0
  
  # R-hat on best chains
  draws_best <- posterior::subset_draws(posterior::as_draws_array(fit),
                                        chain = best_chain_indices)
  rhats_best <- sapply(params_scalar, function(p)
    posterior::rhat(posterior::subset_draws(draws_best, variable = p)))
  
  # Extract filtered samples
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
  
  posterior_samples <- posterior_samples_best
  use_best_chains   <- TRUE
  
}

################################################################################
# 5. PARAMETER ESTIMATES
################################################################################

labels_rho <- c("  0–5  yrs", "  5–18 yrs", "  >18  yrs")

################################################################################
# 6. LOO-CV
################################################################################

log_lik_mat <- if (use_best_chains) {
  posterior_samples$log_lik
} else {
  loo::extract_log_lik(fit, parameter_name = "log_lik")
}

loo_result <- loo(log_lik_mat, save_psis = TRUE, cores = 4)
print(loo_result)

pareto_k <- loo_result$diagnostics$pareto_k

################################################################################
# 7. SAVE RESULTS
################################################################################

out_dir  <- file.path("results", antigen_name,
                      sprintf("skew_normal_fixed_r%.0f_1mu_rho_age", r_fixed_value))
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(fit,              file.path(out_dir, "fit.rds"))
saveRDS(posterior_samples, file.path(out_dir, "posterior_samples.rds"))
saveRDS(loo_result,        file.path(out_dir, "loo_result.rds"))
saveRDS(stan_data,         file.path(out_dir, "stan_data.rds"))

if (use_best_chains) {
  saveRDS(posterior_samples_best, file.path(out_dir, "posterior_samples_best_chains.rds"))
}

summary_df           <- as.data.frame(summary_table)
summary_df$parameter <- rownames(summary_df)
write.csv(summary_df, file.path(out_dir, "parameter_summary.csv"), row.names = FALSE)

################################################################################
# 8. DIAGNOSTIC PLOTS
################################################################################

# Trace plots
pdf(file.path(plot_dir, "trace_plots.pdf"), width = 14, height = 14)
print(traceplot(fit, pars = params, inc_warmup = FALSE))
dev.off()

# Posterior densities
pdf(file.path(plot_dir, "posterior_densities.pdf"), width = 14, height = 14)
print(stan_dens(fit, pars = params, separate_chains = TRUE))
dev.off()

# Posterior intervals
pdf(file.path(plot_dir, "posterior_intervals.pdf"), width = 10, height = 8)
print(mcmc_intervals(as.array(fit, pars = params)))
dev.off()

# Pairs: epidemiological
pdf(file.path(plot_dir, "pairs_epi.pdf"), width = 12, height = 12)
pairs(fit, pars = c("lambda", "alpha", params_rho))
dev.off()

# Pairs: MFI
pdf(file.path(plot_dir, "pairs_mfi.pdf"), width = 12, height = 12)
pairs(fit, pars = c("mu_neg", "delta_sero", "mu_pos",
                    "sigma_neg", "sigma_pos", "skew_neg", "skew_pos"))
dev.off()

# Pareto k
pdf(file.path(plot_dir, "pareto_k_diagnostic.pdf"), width = 10, height = 6)
plot(loo_result, diagnostic = "k", label_points = TRUE)
dev.off()

################################################################################
# 9. CHAIN LOG-LIKELIHOOD DIAGNOSTIC
################################################################################

model_tag <- "4-comp rho-age"
p_chain   <- plot_chain_loglik(fit, antigen_name, model_tag)
ggsave(file.path(plot_dir, "chain_loglik_comparison.png"), p_chain, width = 10, height = 5, dpi = 150)

invisible(fit)
