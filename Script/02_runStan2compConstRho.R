################################################################################
# RUN 2-COMPARTMENT SEROCATALYTIC MODEL — SKEW-NORMAL MFI — DISTRICT
################################################################################

# Run from the repository root (Pv_sero_pcr_modeling/)

if (!exists("antigen_name")) antigen_name <- "EBP"

suppressPackageStartupMessages({
  library(rstan); library(bayesplot); library(loo)
  library(posterior); library(ggplot2); library(dplyr)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

results_dir <- file.path("results", antigen_name, "2compartments_skewnormal")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1. LOAD DATA
################################################################################

data_file <- sprintf("Data/stan_data_%s.rds", antigen_name)
if (!file.exists(data_file)) data_file <- "Data/stan_data.rds"
orig <- readRDS(data_file)

Ab_min <- min(orig$log_mfi)
Ab_max <- max(orig$log_mfi)

stan_data <- list(
  N      = orig$N,
  age    = orig$age,
  Y      = orig$log_mfi,
  Ab_min = Ab_min,
  Ab_max = Ab_max
)

saveRDS(stan_data, file.path(results_dir, "stan_data.rds"))

################################################################################
# 2. COMPILE AND RUN
################################################################################

stan_file <- "Stan/2comp_const_rho.stan"
stan_model_compiled <- stan_model(file = stan_file)

init_fn <- function() list(
  lambda    = runif(1, 0.05, 0.5),
  rho       = runif(1, 0.01, 0.5),
  mu_0      = runif(1, Ab_min + 0.1, Ab_min + (Ab_max - Ab_min) * 0.4),
  mu_1      = runif(1, log(2) + 0.1, (Ab_max - Ab_min) * 0.5),
  epsilon_0 = runif(1, 0.3, 2.0),
  epsilon_1 = runif(1, 0.3, 2.0),
  alpha_0   = runif(1, -2, 5),
  alpha_1   = runif(1, -10, 2)
)

fit <- sampling(
  stan_model_compiled,
  data    = stan_data,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 2,
  seed    = 42,
  #init    = init_fn,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

saveRDS(fit, file.path(results_dir, "fit.rds"))

################################################################################
# 3. CONVERGENCE
################################################################################

params        <- c("lambda", "rho", "mu_0", "mu_1", "epsilon_0", "epsilon_1",
                   "alpha_0", "alpha_1")
summary_table <- summary(fit, pars = params)$summary
rhats         <- summary_table[, "Rhat"]
ess           <- summary_table[, "n_eff"]

print(round(summary_table[, c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")], 4))

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
divergences    <- sum(sapply(sampler_params, function(x) sum(x[,"divergent__"])))

################################################################################
# 4. BEST-CHAINS SELECTION
################################################################################

posterior_samples  <- rstan::extract(fit)
best_chain_indices <- seq_len(fit@sim$chains)
use_best_chains    <- FALSE

if (max(rhats, na.rm = TRUE) > 1.05) {
  ll_arr        <- rstan::extract(fit, pars = "log_lik", permuted = FALSE)
  n_chains      <- dim(ll_arr)[2]
  chain_mean_ll <- sapply(seq_len(n_chains), function(ch) mean(rowSums(ll_arr[, ch, ])))
  
  sorted_order <- order(chain_mean_ll)
  sorted_ll    <- chain_mean_ll[sorted_order]
  gaps         <- diff(sorted_ll)
  gap_pos      <- which.max(gaps)
  gap_size     <- gaps[gap_pos]
  
  best_chain_indices <- sort(sorted_order[(gap_pos + 1):n_chains])
  dropped_chains     <- sort(sorted_order[seq_len(gap_pos)])
  
  within_range <- if (length(best_chain_indices) > 1)
    diff(range(chain_mean_ll[best_chain_indices])) else 0
  ll_diff      <- mean(chain_mean_ll[best_chain_indices]) -
    mean(chain_mean_ll[dropped_chains])
  
  params_filter <- c(params, "prob_seropositive")
  psc <- list()
  for (p in params_filter) {
    arr <- rstan::extract(fit, pars = p, permuted = FALSE)
    sub <- arr[, best_chain_indices, , drop = FALSE]; d <- dim(sub)
    psc[[p]] <- if (length(d)==3 && d[3]==1) as.vector(sub)
    else matrix(sub, nrow=d[1]*d[2], ncol=d[3])
  }
  ll_sub <- ll_arr[, best_chain_indices, , drop = FALSE]
  psc[["log_lik"]] <- matrix(ll_sub, nrow=dim(ll_sub)[1]*dim(ll_sub)[2], ncol=dim(ll_sub)[3])
  saveRDS(psc, file.path(results_dir, "posterior_samples_best_chains.rds"))
  saveRDS(posterior::subset_draws(posterior::as_draws_array(fit), chain = best_chain_indices),
          file.path(results_dir, "fit_best_chains.rds"))
  posterior_samples <- psc
  use_best_chains   <- TRUE
}

saveRDS(posterior_samples, file.path(results_dir, "posterior_samples.rds"))

################################################################################
# 3b. RECOMPUTE CONVERGENCE USING BEST CHAINS ONLY
################################################################################

draws_bc <- posterior::subset_draws(
  posterior::as_draws_array(fit),
  chain    = best_chain_indices,
  variable = params
)
sum_bc <- posterior::summarise_draws(
  draws_bc,
  mean     = mean,
  sd       = sd,
  `2.5%`   = ~quantile(.x, 0.025),
  `25%`    = ~quantile(.x, 0.25),
  `50%`    = ~quantile(.x, 0.50),
  `75%`    = ~quantile(.x, 0.75),
  `97.5%`  = ~quantile(.x, 0.975),
  rhat     = posterior::rhat,
  ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail
)
rhats_bc <- setNames(sum_bc$rhat,     sum_bc$variable)
ess_bc   <- setNames(sum_bc$ess_bulk, sum_bc$variable)

print(as.data.frame(sum_bc[, c("variable","mean","sd","2.5%","50%","97.5%","ess_bulk","rhat")]),
      digits = 4, row.names = FALSE)

################################################################################
# 5. PARAMETER SUMMARY CSV
################################################################################

draws_bc_mupos <- posterior::mutate_variables(
  posterior::subset_draws(draws_bc, variable = c("mu_0", "mu_1")),
  mu_pos = mu_0 + mu_1
)
mu_pos_row <- as.data.frame(posterior::summarise_draws(
  posterior::subset_draws(draws_bc_mupos, variable = "mu_pos"),
  mean     = mean,
  sd       = sd,
  `2.5%`   = ~quantile(.x, 0.025),
  `25%`    = ~quantile(.x, 0.25),
  `50%`    = ~quantile(.x, 0.50),
  `75%`    = ~quantile(.x, 0.75),
  `97.5%`  = ~quantile(.x, 0.975),
  rhat     = posterior::rhat,
  ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail
))
param_summary <- rbind(as.data.frame(sum_bc), mu_pos_row)
names(param_summary)[names(param_summary) == "variable"] <- "parameter"
write.csv(param_summary, file.path(results_dir, "parameter_summary.csv"), row.names = FALSE)

################################################################################
# 6. LOO
################################################################################

log_lik_mat <- posterior_samples[["log_lik"]]
loo_result  <- loo(log_lik_mat, save_psis = TRUE, cores = 4)
saveRDS(loo_result, file.path(results_dir, "loo_result.rds"))
print(loo_result)

################################################################################
# 7. DIAGNOSTIC PLOTS
################################################################################

plots_dir <- file.path(results_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Trace plots
pdf(file.path(plots_dir, "trace_plots.pdf"), width = 12, height = 10)
print(traceplot(fit, pars = params, inc_warmup = FALSE))
dev.off()

# 2. Posterior densities (per chain)
pdf(file.path(plots_dir, "posterior_densities.pdf"), width = 12, height = 10)
print(stan_dens(fit, pars = params, separate_chains = TRUE))
dev.off()

# 3. Posterior intervals
pdf(file.path(plots_dir, "posterior_intervals.pdf"), width = 10, height = 8)
print(mcmc_intervals(as.array(fit, pars = params)))
dev.off()

# 4. Pairs plot — epidemiological params
pdf(file.path(plots_dir, "pairs_plot.pdf"), width = 10, height = 10)
pairs(fit, pars = c("lambda", "rho"))
dev.off()

# 5. Pairs plot — MFI params
pdf(file.path(plots_dir, "pairs_plot_mfi.pdf"), width = 12, height = 12)
pairs(fit, pars = c("mu_0", "mu_1", "epsilon_0", "epsilon_1", "alpha_0", "alpha_1"))
dev.off()

# 6. Pareto k diagnostic
pdf(file.path(plots_dir, "pareto_k_diagnostic.pdf"), width = 10, height = 6)
plot(loo_result, diagnostic = "k", label_points = TRUE)
dev.off()

# 7. Chain log-likelihood comparison
ll_arr_all  <- rstan::extract(fit, pars = "log_lik", permuted = FALSE)
n_iter_pl   <- dim(ll_arr_all)[1]
n_chains_pl <- dim(ll_arr_all)[2]

df_ll <- data.frame(
  iteration = rep(seq_len(n_iter_pl), n_chains_pl),
  chain     = factor(rep(seq_len(n_chains_pl), each = n_iter_pl),
                     labels = paste0("Chain ", seq_len(n_chains_pl))),
  total_ll  = as.vector(apply(ll_arr_all, c(1, 2), sum))
)

chain_colours <- c("Chain 1" = "#2166ac", "Chain 2" = "#d6604d",
                   "Chain 3" = "#f4a582", "Chain 4" = "#4393c3")

p_chain <- ggplot(df_ll, aes(x = iteration, y = total_ll, colour = chain)) +
  geom_line(alpha = 0.8, linewidth = 0.5) +
  scale_colour_manual(values = chain_colours) +
  labs(title  = sprintf("Total log-likelihood per chain — District — %s (2-compartment)", antigen_name),
       x      = "Iteration", y = "Total log-likelihood", colour = "Chain") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(plots_dir, "chain_loglik_comparison.png"), p_chain, width = 10, height = 5, dpi = 150)

# Summary
lambda_v <- posterior_samples$lambda
rho_v    <- posterior_samples$rho
pi_eq    <- mean(lambda_v / (lambda_v + rho_v))
