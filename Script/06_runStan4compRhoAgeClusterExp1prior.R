################################################################################
# RUN STAN — 4-COMPARTMENT MODEL — CLUSTER STRATIFICATION — EXP(1) PRIOR
#   Fixed r | Skew-normal MFI | Age-varying rho | 1 mu_pos
#   One lambda per cluster | All other parameters shared
#
# PRIOR-SENSITIVITY FIX: lambda ~ exponential(1.0) instead of exponential(10.0).
# Model: stan_models/serocatalytic_pcr_4compartments_FIXED_R_skewnormal_1mu_rho_age_cluster_exp1prior.stan
# Data:  Data/stan_data_{antigen}_stratified.rds
#
# Age groups: 0–5 (rho[1]) | 5–18 (rho[2]) | >18 (rho[3])
#
# Output: results/{antigen}/skew_normal_fixed_r{r}_1mu_rho_age_cluster_exp1prior/
################################################################################

rm(list = ls())
# Run from the repository root (Pv_sero_pcr_modeling/)

if (!exists("antigen_name"))  antigen_name  <- "RBP2b.P87"
if (!exists("r_fixed_value")) r_fixed_value <- 8
r_est = r_fixed_value

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
  library(sn)
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

data_file <- sprintf("Data/stan_data_%s_stratified.rds", antigen_name)
if (!file.exists(data_file)) stop("Data not found: ", data_file)

stan_data_orig <- readRDS(data_file)

# Exclude clusters with insufficient sample size (clusters 5 and 32)
clusters_to_exclude <- c(5, 32)
exclude_names <- stan_data_orig$cluster_names[clusters_to_exclude]
keep <- !(stan_data_orig$cluster %in% clusters_to_exclude)

# Remap cluster indices to 1..N_cluster_new
old_clusters          <- stan_data_orig$cluster[keep]
remaining_cluster_ids <- sort(unique(old_clusters))
new_cluster_map       <- setNames(seq_along(remaining_cluster_ids), remaining_cluster_ids)
new_clusters          <- as.integer(new_cluster_map[as.character(old_clusters)])

stan_data <- list(
  N          = sum(keep),
  age        = stan_data_orig$age[keep],
  log_mfi    = stan_data_orig$log_mfi[keep],
  pcr        = stan_data_orig$pcr[keep],
  N_cluster  = length(remaining_cluster_ids),
  cluster    = new_clusters,
  r_fixed    = r_fixed_value,
  lambda_max = stan_data_orig$lambda_max
)

cluster_names <- stan_data_orig$cluster_names[remaining_cluster_ids]

cluster_tab <- table(stan_data$cluster)

################################################################################
# 3. COMPILE AND RUN
################################################################################

stan_file <- "Stan/4comp_rho_age_cluster_exp1prior.stan"

start_compile   <- Sys.time()
stan_model_comp <- stan_model(file = stan_file)

log_mfi_min   <- min(stan_data$log_mfi)
log_mfi_max   <- max(stan_data$log_mfi)

# init_fn <- function() {
#   list(
#     lambda     = rep(0.3, stan_data$N_cluster),
#     alpha      = 0.5,
#     rho        = c(0.5, 0.2, 0.1),
#     mu_neg     = 2.0,
#     delta_sero = 5.0,
#     sigma_neg  = 1.0,
#     sigma_pos  = 1.5,
#     skew_neg   = 1.5,
#     skew_pos   = -5.0
#   )
# }

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

params_lambda <- paste0("lambda[", seq_len(stan_data$N_cluster), "]")
params_rho    <- paste0("rho[", 1:3, "]")
params_shared <- c("alpha",
                   "mu_neg", "mu_pos", "delta_sero",
                   "sigma_neg", "sigma_pos", "skew_neg", "skew_pos")
params_all    <- c(params_lambda, params_shared, params_rho)

summary_shared <- summary(fit, pars = c(params_shared, params_rho))$summary
print(round(summary_shared, 4))

summary_lambda <- summary(fit, pars = params_lambda)$summary
df_lambda_smry <- data.frame(
  cluster = cluster_names,
  median  = round(summary_lambda[, "50%"],    4),
  q025    = round(summary_lambda[, "2.5%"],   4),
  q975    = round(summary_lambda[, "97.5%"],  4),
  mean    = round(summary_lambda[, "mean"],   4),
  sd      = round(summary_lambda[, "sd"],     4),
  Rhat    = round(summary_lambda[, "Rhat"],   4),
  n_eff   = round(summary_lambda[, "n_eff"],  0)
)
print(df_lambda_smry)

summary_all <- summary(fit, pars = params_all)$summary
rhats       <- summary_all[, "Rhat"]
ess_bulk    <- summary_all[, "n_eff"]

if (max(rhats, na.rm = TRUE) >= 1.05) {
  print(names(rhats)[rhats > 1.05])
}

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
divergences    <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
treedepth_hits <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 15)))

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
  
  params_to_filter <- c(params_lambda, params_shared, "rho",
                        "lambda", "prob_seropositive", "prob_pcr_positive", "log_lik")
  posterior_samples_best <- list()
  
  for (param in params_to_filter) {
    arr <- tryCatch(
      rstan::extract(fit, pars = param, permuted = FALSE),
      error = function(e) NULL
    )
    if (is.null(arr)) next
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

labels_rho <- c("0-5  yrs", "5-18 yrs", ">18  yrs")
lambda_mat <- posterior_samples$lambda

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
                      sprintf("skew_normal_fixed_r%.0f_1mu_rho_age_cluster_exp1prior", r_fixed_value))
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

#saveRDS(fit,               file.path(out_dir, "fit.rds"))
saveRDS(posterior_samples, file.path(out_dir, "posterior_samples.rds"))
saveRDS(loo_result,        file.path(out_dir, "loo_result.rds"))
saveRDS(stan_data,         file.path(out_dir, "stan_data.rds"))

if (use_best_chains)
  saveRDS(posterior_samples_best, file.path(out_dir, "posterior_samples_best_chains.rds"))

summary_df           <- as.data.frame(summary_all)
summary_df$parameter <- rownames(summary_df)
write.csv(summary_df, file.path(out_dir, "parameter_summary.csv"), row.names = FALSE)

# Lambda per cluster CSV
write.csv(df_lambda_smry, file.path(out_dir, "lambda_per_cluster.csv"), row.names = FALSE)

################################################################################
# 8. DIAGNOSTIC PLOTS
################################################################################

# Trace: shared parameters
pdf(file.path(plot_dir, "trace_shared.pdf"), width = 14, height = 14)
print(traceplot(fit, pars = c(params_shared, params_rho), inc_warmup = FALSE))
dev.off()

# Trace: lambda per cluster
pdf(file.path(plot_dir, "trace_lambda.pdf"), width = 14, height = 2 * stan_data$N_cluster)
print(traceplot(fit, pars = params_lambda, inc_warmup = FALSE))
dev.off()

# Posterior intervals: shared
pdf(file.path(plot_dir, "posterior_intervals_shared.pdf"), width = 10, height = 8)
print(mcmc_intervals(as.array(fit, pars = c(params_shared, params_rho))))
dev.off()

# Posterior intervals: lambda per cluster
pdf(file.path(plot_dir, "posterior_intervals_lambda.pdf"), width = 10, height = 10)
print(mcmc_intervals(as.array(fit, pars = params_lambda),
                     point_est = "median"))
dev.off()

# Pairs: epi (first 3 lambda + alpha + rho)
n_show <- min(3, stan_data$N_cluster)
pdf(file.path(plot_dir, "pairs_epi.pdf"), width = 12, height = 12)
pairs(fit, pars = c(params_lambda[seq_len(n_show)], "alpha", params_rho))
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

model_tag <- "4-comp rho-age (cluster, exp1 prior)"

p_chain <- plot_chain_loglik(fit, antigen_name, model_tag)
ggsave(file.path(plot_dir, "chain_loglik_comparison.png"),
       p_chain, width = 10, height = 5, dpi = 150)

invisible(fit)
