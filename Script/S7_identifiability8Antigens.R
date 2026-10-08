################################################################################
# IDENTIFIABILITY — all 8 antigens, each antigen's own real fitted parameters,
# own real age distribution, uniform production MCMC settings
#
# chains=4, iter=2000, warmup=1000, adapt_delta=0.99, max_treedepth=15,
# default init, chain selection by log-likelihood if max Rhat > 1.05
#
# Output: results/identifiability_8_antigens/
################################################################################



suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(dplyr)
  library(sn)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("compartment_model.R")

out_dir <- "results/identifiability_8_antigens"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

r_fixed <- 8

antigens <- list(
  
  EBP = list(lambda = 0.5509, alpha = 0.9465, rho = c(3.1007, 2.6695, 0.4683),
             mu_neg = 2.8562, sigma_neg = 1.6217, skew_neg = 5.3783,
             mu_pos = 10.224, sigma_pos = 2.6054, skew_pos = -18.8142),
  
  LF005 = list(lambda = 0.6265, alpha = 0.9802, rho = c(2.1805, 1.1330, 0.3270),
               mu_neg = 3.7028, sigma_neg = 1.5430, skew_neg = 4.4230,
               mu_pos = 7.7325, sigma_pos = 1.4217, skew_pos = -0.9756),
  
  LF010 = list(lambda = 0.5425, alpha = 0.8899, rho = c(3.1969, 2.6299, 0.3960),
               mu_neg = 4.9122, sigma_neg = 0.5643, skew_neg = 5.7920,
               mu_pos = 8.9251, sigma_pos = 1.8440, skew_pos = -1.6765),
  
  LF016 = list(lambda = 0.6314, alpha = 0.8821, rho = c(1.8808, 1.3118, 0.6517),
               mu_neg = 4.6468, sigma_neg = 1.5740, skew_neg = 4.2008,
               mu_pos = 10.2282, sigma_pos = 2.3654, skew_pos = -10.5830),
  
  MSP8 = list(lambda = 0.6258, alpha = 0.9676, rho = c(1.9701, 1.2974, 0.1651),
              mu_neg = 4.5027, sigma_neg = 1.5327, skew_neg = 9.7692,
              mu_pos = 8.5303, sigma_pos = 1.4324, skew_pos = -1.1687),
  
  PTEX150 = list(lambda = 0.6902, alpha = 0.0607, rho = c(0.6188, 0.0971, 0.0246),
                 mu_neg = 4.4827, sigma_neg = 1.5851, skew_neg = 4.6255,
                 mu_pos = 8.4931, sigma_pos = 1.6673, skew_pos = -2.3823),
  
  PvCSS = list(lambda = 0.601, alpha = 0.959, rho = c(2.427, 1.729, 0.361),
               mu_neg = 2.982, sigma_neg = 1.248, skew_neg = 5.526,
               mu_pos = 7.250, sigma_pos = 1.488, skew_pos = -0.399),
  
  RBP2b = list(lambda = 0.5711, alpha = 0.9619, rho = c(2.7209, 2.2892, 0.7725),
               mu_neg = 4.4710, sigma_neg = 0.9292, skew_neg = 4.0095,
               mu_pos = 10.2048, sigma_pos = 1.5125, skew_pos = -19.5527)
)

antigen_file_names <- c(
  EBP = "EBP", LF005 = "LF005", LF010 = "LF010", LF016 = "LF016",
  MSP8 = "MSP8", PTEX150 = "PTEX150", PvCSS = "PvCSS", RBP2b = "RBP2b.P87"
)

stan_model_comp <- stan_model(file = "stan_models/serocatalytic_pcr_4compartments_FIXED_R_skewnormal_1mu_rho_age.stan")

all_recovery <- list()

for (ag in names(antigens)) {
  
  p <- antigens[[ag]]
  
  true_vals <- c(
    lambda = p$lambda, alpha = p$alpha,
    `rho[1]` = p$rho[1], `rho[2]` = p$rho[2], `rho[3]` = p$rho[3],
    mu_neg = p$mu_neg, mu_pos = p$mu_pos,
    sigma_neg = p$sigma_neg, sigma_pos = p$sigma_pos,
    skew_neg = p$skew_neg, skew_pos = p$skew_pos
  )
  
  sd_real <- readRDS(sprintf("data/stan_data_%s.rds", antigen_file_names[ag]))
  ages <- sd_real$age
  N <- length(ages)
  lambda_max <- sd_real$lambda_max
  
  set.seed(2024)
  pcr_vec <- integer(N)
  log_mfi <- numeric(N)
  
  for (i in seq_len(N)) {
    S <- calculate_compartments_age_varying(ages[i], p$lambda, p$alpha, p$rho, r_fixed)
    state <- sample(1:4, 1, prob = pmax(S, 0))
    pcr_vec[i] <- as.integer(state %in% c(2, 3))
    if (state %in% c(1, 2)) {
      log_mfi[i] <- rsn(1, xi = p$mu_neg, omega = p$sigma_neg, alpha = p$skew_neg)
    } else {
      log_mfi[i] <- rsn(1, xi = p$mu_pos, omega = p$sigma_pos, alpha = p$skew_pos)
    }
  }
  
  
  
  stan_data <- list(
    N = N, age = ages, log_mfi = log_mfi, pcr = pcr_vec,
    r_fixed = r_fixed, lambda_max = lambda_max
  )
  
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
  rhats <- summary_table[, "Rhat"]
  
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  div <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  treedepth_hits <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 15)))
  
  
  posterior_samples <- rstan::extract(fit)
  use_best_chains <- FALSE
  best_chain_indices <- seq_len(fit@sim$chains)
  chains_dropped <- 0
  
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
    chains_dropped <- length(dropped_chains)
    ll_diff <- mean(chain_mean_ll[best_chain_indices]) -
      mean(chain_mean_ll[dropped_chains])
    
    within_range <- if (length(best_chain_indices) > 1)
      diff(range(chain_mean_ll[best_chain_indices])) else 0
    
    
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
    
    posterior_samples <- posterior_samples_best
    use_best_chains <- TRUE
    cat(">>> Continuing with BEST CHAINS ONLY <<<\n\n")
    
  } else {
    cat("All chains converged well. Using all chains.\n\n")
  }
  
  ps <- posterior_samples
  pooled <- list(
    lambda = ps$lambda, alpha = ps$alpha,
    `rho[1]` = ps$rho[,1], `rho[2]` = ps$rho[,2], `rho[3]` = ps$rho[,3],
    mu_neg = ps$mu_neg, mu_pos = ps$mu_pos,
    sigma_neg = ps$sigma_neg, sigma_pos = ps$sigma_pos,
    skew_neg = ps$skew_neg, skew_pos = ps$skew_pos
  )
  
  recovery_df <- data.frame(antigen = ag,
                            parameter = names(true_vals), true_value = unname(true_vals),
                            stringsAsFactors = FALSE)
  recovery_df$median <- NA; recovery_df$q2.5 <- NA; recovery_df$q97.5 <- NA
  
  for (i in seq_len(nrow(recovery_df))) {
    samps <- pooled[[recovery_df$parameter[i]]]
    recovery_df$median[i] <- median(samps)
    recovery_df$q2.5[i]   <- quantile(samps, 0.025)
    recovery_df$q97.5[i]  <- quantile(samps, 0.975)
  }
  
  recovery_df$recovered <- ifelse(
    recovery_df$true_value >= recovery_df$q2.5 & recovery_df$true_value <= recovery_df$q97.5,
    "YES", "NO"
  )
  recovery_df$rel_bias <- round((recovery_df$median - recovery_df$true_value) / recovery_df$true_value * 100, 1)
  recovery_df$max_rhat <- max(rhats, na.rm = TRUE)
  recovery_df$divergences <- div
  recovery_df$treedepth_hits <- treedepth_hits
  recovery_df$chains_dropped <- chains_dropped
  recovery_df$use_best_chains <- use_best_chains
  
  print(recovery_df[, c("parameter","true_value","median","q2.5","q97.5","recovered","rel_bias")])
  cat("\n")
  
  write.csv(recovery_df, file.path(out_dir, sprintf("recovery_%s.csv", ag)), row.names = FALSE)
  all_recovery[[ag]] <- recovery_df
  
  rm(fit); gc()
}

combined <- bind_rows(all_recovery)
write.csv(combined, file.path(out_dir, "recovery_8_antigens.csv"), row.names = FALSE)

