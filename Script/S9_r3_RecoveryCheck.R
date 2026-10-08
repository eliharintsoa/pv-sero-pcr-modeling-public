################################################################################
# PARAMETER RECOVERY AT r = 3/year — EBP and LF016
#
# True parameters = each antigen's own fitted values AT r=3 (EBP: from the
# existing r-sensitivity grid; LF016: from S11_lf016_real_fit_r3.R), real
# ages, real N, production MCMC settings. Tests whether r=3 recovers lambda
# (and alpha, rho) better than r=8 did for EBP (which failed at r=8) and
# whether it still recovers for LF016 (which passed at r=8).
#
# Output: results/identifiability_r3/
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

out_dir <- "results/identifiability_r3"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

r_fixed <- 3

lf016_fit <- read.csv("results/LF016/skew_normal_fixed_r3_1mu_rho_age/fitted_medians.csv")
lf016_vals <- setNames(lf016_fit$value, lf016_fit$parameter)

antigens <- list(
  EBP = list(lambda = 0.230933, alpha = 0.956321, rho = c(1.907924, 1.345370, 0.192977),
             mu_neg = 2.8562, sigma_neg = 1.6217, skew_neg = 5.3783,
             mu_pos = 10.224, sigma_pos = 2.6054, skew_pos = -18.8142),
  
  LF016 = list(lambda = unname(lf016_vals["lambda"]), alpha = unname(lf016_vals["alpha"]),
               rho = unname(c(lf016_vals["rho1"], lf016_vals["rho2"], lf016_vals["rho3"])),
               mu_neg = unname(lf016_vals["mu_neg"]), sigma_neg = unname(lf016_vals["sigma_neg"]),
               skew_neg = unname(lf016_vals["skew_neg"]), mu_pos = unname(lf016_vals["mu_pos"]),
               sigma_pos = unname(lf016_vals["sigma_pos"]), skew_pos = unname(lf016_vals["skew_pos"]))
)

antigen_file_names <- c(EBP = "EBP",
                        LF016 = "LF016")

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
    
    sorted_order <- order(chain_mean_ll)
    gaps <- diff(chain_mean_ll[sorted_order])
    gap_pos <- which.max(gaps)
    best_chain_indices <- sort(sorted_order[(gap_pos + 1):n_chains_total])
    dropped_chains <- sort(sorted_order[seq_len(gap_pos)])
    chains_dropped <- length(dropped_chains)
    
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
  recovery_df$chains_dropped <- chains_dropped
  
  print(recovery_df[, c("parameter","true_value","median","q2.5","q97.5","recovered","rel_bias")])
  cat("\n")
  
  write.csv(recovery_df, file.path(out_dir, sprintf("recovery_%s_r3.csv", ag)), row.names = FALSE)
  all_recovery[[ag]] <- recovery_df
  
  rm(fit); gc()
}

combined <- bind_rows(all_recovery)
write.csv(combined, file.path(out_dir, "recovery_r3_combined.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s/recovery_r3_combined.csv\n", out_dir))

# Plot: lambda/alpha/rho recovery, r=8 vs r=3, EBP vs LF016
library(ggplot2)

plot_files <- list(
  list(antigen = "EBP",   r = "r=8", path = "results/identifiability_8_antigens/recovery_EBP.csv"),
  list(antigen = "EBP",   r = "r=3", path = file.path(out_dir, "recovery_EBP_r3.csv")),
  list(antigen = "LF016", r = "r=8", path = "results/identifiability_8_antigens/recovery_LF016.csv"),
  list(antigen = "LF016", r = "r=3", path = file.path(out_dir, "recovery_LF016_r3.csv"))
)

keep_params <- c("lambda", "alpha", "rho[1]", "rho[2]", "rho[3]")

plot_df <- bind_rows(lapply(plot_files, function(f) {
  d <- read.csv(f$path)
  d <- d[d$parameter %in% keep_params, ]
  d$antigen_label <- f$antigen
  d$r_label <- f$r
  d
}))

plot_df$parameter <- factor(plot_df$parameter, levels = keep_params)
plot_df$r_label <- factor(plot_df$r_label, levels = c("r=8", "r=3"))
plot_df$antigen_label <- factor(plot_df$antigen_label, levels = c("EBP", "LF016"),
                                labels = c("EBP", "MSP1-19"))

p <- ggplot(plot_df, aes(x = r_label, y = median, colour = antigen_label)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.15,
                position = position_dodge(width = 0.4), linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.4), size = 2.5) +
  geom_point(aes(y = true_value, shape = "True value"),
             position = position_dodge(width = 0.4), size = 2.8) +
  facet_wrap(~ parameter, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c("EBP" = "#D55E00", "MSP1-19" = "#0072B2"), name = "Antigen") +
  scale_shape_manual(values = c("True value" = 18), name = NULL) +
  labs(x = NULL, y = "Posterior median [95% CrI]") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

dir.create("paper", showWarnings = FALSE)
ggsave("paper/figure_r3_vs_r8_recovery.png", p, width = 11, height = 4, dpi = 300, bg = "white")
ggsave("paper/figure_r3_vs_r8_recovery.pdf", p, width = 11, height = 4)

