
################################################################################
# RUN STAN MODEL WITH FIXED R - SKEW-NORMAL MFI DISTRIBUTION
# Uses prepared data and fixes r from literature
################################################################################

# Clear workspace
rm(list = ls())

# Set working directory
# Run from the repository root (Pv_sero_pcr_modeling/)

# Antigen name (can be overridden before sourcing)

if (!exists("antigen_name")) antigen_name <- "MSP8"


# r_fixed: PCR clearance rate from literature (per year)

if (!exists("r_fixed_value")) {
  r_fixed_value <- 8
}

################################################################################
# 1. LOAD LIBRARIES
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

# Configure Stan
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
# 2. LOAD PREPARED DATA
################################################################################

data_file <- sprintf("Data/stan_data_%s.rds", antigen_name)
if (!file.exists(data_file)) {
  data_file <- "Data/stan_data.rds"
}

if (!file.exists(data_file)) {
  stop("Data preparation required")
}

stan_data_orig <- readRDS(data_file)

# Build stan_data for fixed-r model:
# Replace r_max with r_fixed, keep everything else
stan_data <- list(
  N        = stan_data_orig$N,
  age      = stan_data_orig$age,
  log_mfi  = stan_data_orig$log_mfi,
  pcr      = stan_data_orig$pcr,
  r_fixed  = r_fixed_value,
  lambda_max = stan_data_orig$lambda_max
)

################################################################################
# 3. COMPILE AND RUN STAN MODEL
################################################################################

stan_file <- "Stan/4comp_const_rho.stan"

# Compilation
start_compile <- Sys.time()

stan_model_compiled <- stan_model(file = stan_file)

end_compile <- Sys.time()

# Init function (no r parameter)
log_mfi_min <- min(stan_data$log_mfi)
log_mfi_max <- max(stan_data$log_mfi)
log_mfi_range <- log_mfi_max - log_mfi_min

init_fn <- function() {
  mu_neg_init <- runif(1, max(log_mfi_min + 0.1, 2.5), min(5.0, log_mfi_max - 4.5))
  delta_mu_init <- runif(1, max(4.1, 4.5), min(log_mfi_range - 0.1, 7.0))
  list(
    lambda    = runif(1, 0.05, 0.3),
    alpha     = runif(1, 0.3, 0.9),
    rho       = runif(1, 0.01, 0.5),
    mu_neg    = mu_neg_init,
    delta_mu  = delta_mu_init,
    sigma_neg = runif(1, 0.5, 2.0),
    sigma_pos = runif(1, 0.5, 2.0),
    skew_neg  = runif(1, -2, 5),
    skew_pos  = runif(1, -10, 2)
  )
}

# Sampling
start_sampling <- Sys.time()

fit <- sampling(
  stan_model_compiled,
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  thin = 1,
  cores = 4,
  seed = 42,
  #init = init_fn,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

end_sampling <- Sys.time()
sampling_time <- as.numeric(end_sampling - start_sampling, units = "mins")

################################################################################
# 4. CONVERGENCE DIAGNOSTICS
################################################################################

# Parameters of interest (NO r - it's fixed)
params <- c("lambda", "alpha", "rho",
            "mu_neg", "mu_pos", "sigma_neg", "sigma_pos",
            "skew_neg", "skew_pos")

# Summary statistics
summary_table <- summary(fit, pars = params)$summary
print(summary_table)

# Check R-hat
rhats <- summary_table[, "Rhat"]
if (max(rhats, na.rm = TRUE) >= 1.05) {
  print(names(rhats)[rhats > 1.05])
}

# Check Effective Sample Size
ess_bulk <- summary_table[, "n_eff"]
if (min(ess_bulk, na.rm = TRUE) <= 400) {
  print(names(ess_bulk)[ess_bulk < 400])
}

# Divergent transitions
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))

# Tree depth
max_treedepth_hits <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 15)))

################################################################################
# 4b. CHAIN SELECTION: KEEP BEST CHAINS IF CONVERGENCE IS POOR
################################################################################

# Extract all posterior samples first (used below and potentially replaced)
posterior_samples  <- rstan::extract(fit)
use_best_chains    <- FALSE
best_chain_indices <- seq_len(fit@sim$chains)

if (max(rhats, na.rm = TRUE) > 1.05) {
  
  ll_arr         <- rstan::extract(fit, pars = "log_lik", permuted = FALSE)
  n_chains_total <- dim(ll_arr)[2]
  
  chain_mean_ll <- sapply(seq_len(n_chains_total), function(ch) {
    mean(rowSums(ll_arr[, ch, ]))
  })
  
  # --- Find the natural split by detecting the largest gap in sorted ll ---
  # Sort chains by mean log-likelihood (ascending)
  sorted_order <- order(chain_mean_ll)
  sorted_ll    <- chain_mean_ll[sorted_order]
  
  # Gaps between consecutive sorted values
  gaps     <- diff(sorted_ll)
  gap_pos  <- which.max(gaps)   # position of largest gap in sorted sequence
  gap_size <- gaps[gap_pos]     # magnitude of that gap
  
  # Chains above the gap -> better mode (keep). below -> worse (drop)
  best_chain_indices <- sort(sorted_order[(gap_pos + 1):n_chains_total])
  dropped_chains     <- sort(sorted_order[seq_len(gap_pos)])
  
  # Within-group spread of kept chains (should be small if they mixed well)
  within_range <- if (length(best_chain_indices) > 1)
    diff(range(chain_mean_ll[best_chain_indices])) else 0
  ll_diff      <- mean(chain_mean_ll[best_chain_indices]) -
    mean(chain_mean_ll[dropped_chains])
  
  # Extract filtered posterior samples (all params + prob_ quantities)
  params_to_filter <- c(params, "prob_seropositive", "prob_pcr_positive")
  posterior_samples_best_chains <- list()
  
  for (param in params_to_filter) {
    arr <- rstan::extract(fit, pars = param, permuted = FALSE)
    sub <- arr[, best_chain_indices, , drop = FALSE]
    d   <- dim(sub)
    if (length(d) == 3 && d[3] == 1) {
      posterior_samples_best_chains[[param]] <- as.vector(sub)
    } else {
      posterior_samples_best_chains[[param]] <- matrix(sub,
                                                       nrow = d[1] * d[2], ncol = d[3])
    }
  }
  
  # Carry log_lik in filtered samples (needed for LOO)
  ll_sub <- ll_arr[, best_chain_indices, , drop = FALSE]
  posterior_samples_best_chains[["log_lik"]] <- matrix(ll_sub,
                                                       nrow = dim(ll_sub)[1] * dim(ll_sub)[2], ncol = dim(ll_sub)[3])
  
  # Check R-hat on best chains alone
  draws_best <- posterior::subset_draws(posterior::as_draws_array(fit),
                                        chain = best_chain_indices)
  rhats_best <- sapply(params, function(p) {
    posterior::rhat(posterior::subset_draws(draws_best, variable = p))
  })
  
  # Replace posterior_samples for all subsequent analysis
  posterior_samples <- posterior_samples_best_chains
  use_best_chains   <- TRUE
  
}

################################################################################
# 5. PARAMETER ESTIMATES
################################################################################

lambda_mean <- mean(posterior_samples$lambda)
alpha_mean <- mean(posterior_samples$alpha)
rho_mean <- mean(posterior_samples$rho)

################################################################################
# 6. MODEL COMPARISON (LOO-CV)
################################################################################

if (use_best_chains) {
  log_lik <- posterior_samples$log_lik        # already filtered to best chains
} else {
  log_lik <- extract_log_lik(fit, parameter_name = "log_lik")
}

loo_result <- loo(log_lik, save_psis = TRUE, cores = 4)
print(loo_result)

pareto_k <- loo_result$diagnostics$pareto_k
n_high_k <- sum(pareto_k > 0.7)


################################################################################
# 7. SAVE RESULTS
################################################################################

out_dir <- file.path("results", antigen_name, sprintf("skew_normal_fixed_r%.0f", r_fixed_value))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

saveRDS(fit, file.path(out_dir, "fit_fixed_r.rds"))
saveRDS(posterior_samples, file.path(out_dir, "posterior_samples.rds"))

if (use_best_chains) {
  saveRDS(posterior_samples_best_chains,
          file.path(out_dir, "posterior_samples_best_chains.rds"))
  saveRDS(posterior::subset_draws(posterior::as_draws_array(fit),
                                  chain = best_chain_indices),
          file.path(out_dir, "fit_best_chains.rds"))
}

saveRDS(loo_result, file.path(out_dir, "loo_result.rds"))
saveRDS(stan_data, file.path(out_dir, "stan_data.rds"))

summary_df <- as.data.frame(summary_table)
summary_df$parameter <- rownames(summary_df)
write.csv(summary_df, file.path(out_dir, "parameter_summary.csv"), row.names = FALSE)

################################################################################
# 8. DIAGNOSTIC PLOTS
################################################################################

plot_dir <- file.path(out_dir, "plots")

# 1. Trace plots
pdf(file.path(plot_dir, "trace_plots.pdf"), width = 12, height = 12)
print(traceplot(fit, pars = params, inc_warmup = FALSE))
dev.off()

# 2. Density plots
pdf(file.path(plot_dir, "posterior_densities.pdf"), width = 12, height = 12)
print(stan_dens(fit, pars = params, separate_chains = TRUE))
dev.off()

# 3. Posterior intervals
pdf(file.path(plot_dir, "posterior_intervals.pdf"), width = 10, height = 8)
print(mcmc_intervals(as.array(fit, pars = params)))
dev.off()

# 4. Pairs plot (epidemiological params - no r)
pdf(file.path(plot_dir, "pairs_plot.pdf"), width = 12, height = 12)
pairs(fit, pars = c("lambda", "alpha", "rho"))
dev.off()

# 5. Pairs plot (MFI params)
pdf(file.path(plot_dir, "pairs_plot_mfi.pdf"), width = 12, height = 12)
pairs(fit, pars = c("mu_neg", "mu_pos", "sigma_neg", "sigma_pos", "skew_neg", "skew_pos"))
dev.off()

# 6. Pareto k diagnostic
pdf(file.path(plot_dir, "pareto_k_diagnostic.pdf"), width = 10, height = 6)
plot(loo_result, diagnostic = "k", label_points = TRUE)
dev.off()

################################################################################
# 9. CHAIN LOG-LIKELIHOOD DIAGNOSTIC
################################################################################

model_tag <- "4-comp const-rho"
p_chain   <- plot_chain_loglik(fit, antigen_name, model_tag)
ggsave(file.path(plot_dir, "chain_loglik_comparison.png"), p_chain, width = 10, height = 5, dpi = 150)

invisible(fit)
