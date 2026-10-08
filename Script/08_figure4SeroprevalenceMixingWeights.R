################################################################################
# COMPUTE PCR-CONDITIONED MIXING WEIGHTS — ALL FOUR MODELS
# + FIGURE 5 — age-stratified seroprevalence by PCR status (4C-rho-age only)
#
# For each antigen × model × age, computes:
#
#   4C models (PCR data used):
#     w_neg(a) = P0(a) / (N0(a) + P0(a))  = P(seropos | PCR-, age a)
#     w_pos(a) = P1(a) / (N1(a) + P1(a))  = P(seropos | PCR+, age a)
#     pi(a)    = P1(a) + P0(a)             = marginal seroprevalence
#
#   2C models (no PCR data):
#     w_neg = w_pos = pi(a) = seroprevalence at age a
#     (no distinction between PCR groups — single weight for all individuals)
#
# Constant-rho models use best-chain posteriors when available.
#
# The Figure 5 age-stratified seroprevalence plot (0-5 / 6-18 / 18+ yr, 4C-rho-age
# only) is built from the same age-group summary used for the mixing-weights
# table below, so the figure and the table are guaranteed to be the same
# underlying numbers.
#
# Outputs:
#   paper/tables/mixing_weights_summary.csv       median [95% CrI] per model × antigen × age
#   paper/tables/mixing_weights_draws.csv         all draws (long format, for plotting)
#   paper/tables/mixing_weights_agegroup.csv      <18/>=18 yr summary, all models
#   paper/tables/seroprevalence_agegroup_pcr.csv   0-5/6-18/18+ yr summary, 4C-rho-age only
#                                                     (for text: matches the model's own age bands)
#   paper/tables/seroprevalence_finebin_pcr.csv    9-bin summary, 4C-rho-age only
#                                                     (for Figure 5: 0-5/5-10/10-15/15-20/20-30/
#                                                     30-40/40-50/50-60/60+)
#   paper/seroprevalence_observed_vs_predicted_finebin_age         ( Figure 4 )
################################################################################

# Run from the repository root (Pv_sero_pcr_modeling/)
source("compartment_model.R")   # calculate_compartments_age_varying(), a1, a2

suppressPackageStartupMessages({
  library(dplyr)
  library(sn)
  library(ggplot2)
  library(scales)
})

r_fixed  <- 8.0
N_DRAWS  <- 4000
KEY_AGES <- c(2, 5, 10, 15, 18, 25, 30, 50, 70)

antigens <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")
ANTIGEN_COL <- c(EBP = "log_mfi_EBP", LF005 = "log_mfi_LF005", LF010 = "log_mfi_LF010",
                 LF016 = "log_mfi_LF016", MSP8 = "log_mfi_MSP8",
                 PTEX150 = "log_mfi_PTEX150", PvCSS = "log_mfi_PvCSS",
                 RBP2b.P87 = "log_mfi_RBP2b.P87")
ANTIGEN_LABELS <- c(EBP = "EBP", LF005 = "Pv-fam-a", LF010 = "MSP5", LF016 = "MSP1-19",
                    MSP8 = "MSP8", PTEX150 = "PTEX150", PvCSS = "PvCSS",
                    RBP2b.P87 = "RBP2b")

models <- list(
  "2C-const"   = list(dir = "2compartments_skewnormal",
                      n_comp = 2, rho_age = FALSE, best_chains = TRUE),
  "2C-rho-age" = list(dir = "skew_normal_2comp_rho_age",
                      n_comp = 2, rho_age = TRUE,  best_chains = FALSE),
  "4C-const"   = list(dir = "skew_normal_fixed_r8",
                      n_comp = 4, rho_age = FALSE, best_chains = TRUE),
  "4C-rho-age" = list(dir = "skew_normal_fixed_r8_1mu_rho_age",
                      n_comp = 4, rho_age = TRUE,  best_chains = FALSE)
)

# ── 2-compartment seroprevalence functions ────────────────────────────────────

# Constant rho: analytic formula P(a) = λ/(λ+ρ) × (1 - exp(-(λ+ρ)a))
seroprev_2c_const <- function(a, lambda, rho) {
  lambda / (lambda + rho) * (1 - exp(-(lambda + rho) * a))
}

# Age-varying rho: chained across boundaries a1=5, a2=18
seroprev_2c_rho_age <- function(a, lambda, rho_vec) {
  # segment 1: [0, a1) with rho1
  P_eq1 <- lambda / (lambda + rho_vec[1])
  if (a <= a1) {
    return(P_eq1 * (1 - exp(-(lambda + rho_vec[1]) * a)))
  }
  P_a1 <- P_eq1 * (1 - exp(-(lambda + rho_vec[1]) * a1))
  
  # segment 2: [a1, a2) with rho2, starting from P_a1
  P_eq2 <- lambda / (lambda + rho_vec[2])
  if (a <= a2) {
    return(P_eq2 + (P_a1 - P_eq2) * exp(-(lambda + rho_vec[2]) * (a - a1)))
  }
  P_a2 <- P_eq2 + (P_a1 - P_eq2) * exp(-(lambda + rho_vec[2]) * (a2 - a1))
  
  # segment 3: [a2, inf) with rho3, starting from P_a2
  P_eq3 <- lambda / (lambda + rho_vec[3])
  P_eq3 + (P_a2 - P_eq3) * exp(-(lambda + rho_vec[3]) * (a - a2))
}

# ── helpers ───────────────────────────────────────────────────────────────────

ci_pct <- function(x, d = 1) {
  q <- quantile(x * 100, c(0.025, 0.5, 0.975), na.rm = TRUE)
  sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), q[2], q[1], q[3])
}

load_ps <- function(ag, mkey) {
  cfg  <- models[[mkey]]
  dir  <- file.path("results", ag, cfg$dir)
  path <- file.path(dir, "posterior_samples.rds")
  if (cfg$best_chains) {
    bc <- file.path(dir, "posterior_samples_best_chains.rds")
    if (file.exists(bc)) path <- bc
  }
  if (!file.exists(path)) return(NULL)
  ps <- readRDS(path)
  # harmonise rho storage
  if (is.null(dim(ps$rho)) && !is.null(ps[["rho[1]"]]))
    ps$rho <- cbind(ps[["rho[1]"]], ps[["rho[2]"]], ps[["rho[3]"]])
  ps
}

# ── main loop ─────────────────────────────────────────────────────────────────

summary_rows <- list()
draws_rows   <- list()

for (mkey in names(models)) {
  cfg <- models[[mkey]]
  
  for (ag in antigens) {
    
    ps <- load_ps(ag, mkey)
    if (is.null(ps)) next
    
    n   <- min(N_DRAWS, length(ps$lambda))
    lam <- ps$lambda[1:n]
    rho <- if (is.matrix(ps$rho)) ps$rho[1:n, , drop = FALSE] else ps$rho[1:n]
    alp <- if (!is.null(ps$alpha)) ps$alpha[1:n] else rep(NA_real_, n)
    
    # matrices: rows = draws, cols = ages
    w_neg_mat <- matrix(NA_real_, n, length(KEY_AGES))
    w_pos_mat <- matrix(NA_real_, n, length(KEY_AGES))
    pi_mat    <- matrix(NA_real_, n, length(KEY_AGES))
    
    for (i in seq_len(n)) {
      rho_i <- if (is.matrix(rho)) rho[i, ] else rho[i]
      
      for (j in seq_along(KEY_AGES)) {
        a <- KEY_AGES[j]
        
        if (cfg$n_comp == 2) {
          # 2C model: single weight = marginal seroprevalence
          pi_val <- if (cfg$rho_age) {
            seroprev_2c_rho_age(a, lam[i], rho_i)
          } else {
            seroprev_2c_const(a, lam[i], rho_i)
          }
          w_neg_mat[i, j] <- pi_val   # no PCR conditioning: w- = w+ = pi
          w_pos_mat[i, j] <- pi_val
          pi_mat[i, j]    <- pi_val
          
        } else {
          # 4C model: PCR-conditioned weights
          rho_vec_i <- if (length(rho_i) == 1) rep(rho_i, 3) else rho_i
          S  <- calculate_compartments_age_varying(a, lam[i], alp[i], rho_vec_i, r_fixed)
          N0 <- S["N0"]; N1 <- S["N1"]; P1 <- S["P1"]; P0 <- S["P0"]
          w_neg_mat[i, j] <- P0 / max(N0 + P0, 1e-10)
          w_pos_mat[i, j] <- P1 / max(N1 + P1, 1e-10)
          pi_mat[i, j]    <- P1 + P0
        }
      }
    }
    
    # summary
    for (j in seq_along(KEY_AGES)) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        model   = mkey,
        antigen = ag,
        age     = KEY_AGES[j],
        w_neg   = ci_pct(w_neg_mat[, j]),
        w_pos   = ci_pct(w_pos_mat[, j]),
        pi      = ci_pct(pi_mat[, j]),
        pcr_conditioned = cfg$n_comp == 4,
        stringsAsFactors = FALSE
      )
    }
    
    # draws (long format)
    for (j in seq_along(KEY_AGES)) {
      draws_rows[[length(draws_rows) + 1]] <- data.frame(
        model   = mkey,
        antigen = ag,
        age     = KEY_AGES[j],
        draw    = seq_len(n),
        w_neg   = w_neg_mat[, j],
        w_pos   = w_pos_mat[, j],
        pi      = pi_mat[, j],
        pcr_conditioned = cfg$n_comp == 4,
        stringsAsFactors = FALSE
      )
    }
  }
}

# ── save ──────────────────────────────────────────────────────────────────────

dir.create("paper/tables", showWarnings = FALSE, recursive = TRUE)

summary_df <- bind_rows(summary_rows)
write.csv(summary_df, "paper/tables/mixing_weights_summary.csv", row.names = FALSE)

draws_df <- bind_rows(draws_rows)
write.csv(draws_df, "paper/tables/mixing_weights_draws.csv", row.names = FALSE)

# ── age-group summary: 0-5 / 6-18 / 18+ yr ────────────────────────────────────
# For each draw, average weights across KEY_AGES within each group, then
# summarise across draws to get median [95% CrI]. Same three bands used by
#   0-5 yr  -> KEY_AGES {2, 5}
#   6-18 yr -> KEY_AGES {10, 15, 18}
#   18+ yr  -> KEY_AGES {25, 30, 50, 70}

agegroup_draws <- draws_df %>%
  mutate(age_group = case_when(
    age <= 5  ~ "0-5 yr",
    age <= 18 ~ "6-18 yr",
    TRUE      ~ "18+ yr"
  )) %>%
  group_by(model, antigen, pcr_conditioned, age_group, draw) %>%
  summarise(
    w_neg = mean(w_neg, na.rm = TRUE),
    w_pos = mean(w_pos, na.rm = TRUE),
    pi    = mean(pi,    na.rm = TRUE),
    .groups = "drop"
  )

agegroup_summary <- agegroup_draws %>%
  group_by(model, antigen, pcr_conditioned, age_group) %>%
  summarise(
    w_neg = ci_pct(w_neg),
    w_pos = ci_pct(w_pos),
    pi    = ci_pct(pi),
    .groups = "drop"
  ) %>%
  mutate(age_group = factor(age_group, levels = c("0-5 yr", "6-18 yr", "18+ yr"))) %>%
  arrange(antigen, model, age_group)

write.csv(agegroup_summary, "paper/tables/mixing_weights_agegroup.csv", row.names = FALSE)

# ==============================================================================
# FIGURE 5 — age-stratified seroprevalence by PCR status, 4C-rho-age only
#
# ==============================================================================

################################################################################
# OBSERVED (data-driven, Method 3 / p+|x = 0.5) vs PREDICTED seroprevalence,
# https://www.mdpi.com/2227-7390/12/2/217
#
# Output: paper/tables/seroprevalence_observed_vs_predicted_finebin_age.csv
#         paper/seroprevalence_observed_vs_predicted_finebin_age.png
################################################################################



BIN_WIDTH <- 10

d_all <- read.csv("data/data_cleaned.csv")

parse_median <- function(x) as.numeric(sub("^([-0-9.]+).*$", "\\1", x))
params_all <- read.csv("paper/table2_4comp_rho_age_parameters.csv") %>%
  mutate(across(c(mu_neg, mu_pos, sigma_neg, sigma_pos, skew_neg, skew_pos), parse_median))

find_cutoff_method3 <- function(mu_neg, sigma_neg, skew_neg, mu_pos, sigma_pos, skew_pos, pi) {
  f <- function(x) {
    pi * dsn(x, xi = mu_pos, omega = sigma_pos, alpha = skew_pos) -
      (1 - pi) * dsn(x, xi = mu_neg, omega = sigma_neg, alpha = skew_neg)
  }
  lo <- min(mu_neg, mu_pos); hi <- max(mu_neg, mu_pos)
  uniroot(f, interval = c(lo, hi))$root
}

binom_ci <- function(x, n) {
  if (n == 0) return(c(lower = NA_real_, upper = NA_real_))
  ci <- binom.test(x, n)$conf.int
  c(lower = ci[1], upper = ci[2])
}

max_age <- ceiling(max(d_all$agen) / BIN_WIDTH) * BIN_WIDTH
breaks  <- seq(0, max_age, BIN_WIDTH) 

all_obs <- list()
all_pred <- list()

for (ag in antigens) {
  mfi_col <- ANTIGEN_COL[[ag]]
  p <- params_all %>% filter(Antigen == ag)
  
  ps <- readRDS(file.path("results", ag, "skew_normal_fixed_r8_1mu_rho_age", "posterior_samples.rds"))
  stan_data <- readRDS(file.path("results", ag, "skew_normal_fixed_r8_1mu_rho_age", "stan_data.rds"))
  pi_overall <- mean(colMeans(ps$prob_seropositive))
  
  cutoff <- find_cutoff_method3(p$mu_neg, p$sigma_neg, p$skew_neg,
                                p$mu_pos, p$sigma_pos, p$skew_pos, pi_overall)
  
  d <- d_all
  d$seropos  <- as.integer(d[[mfi_col]] > cutoff)
  d$age_bin  <- cut(d$agen, breaks = breaks, include.lowest = TRUE)
  d$age_mid  <- (breaks[as.integer(d$age_bin)] + breaks[as.integer(d$age_bin) + 1]) / 2
  
  tab <- d %>%
    group_by(age_mid) %>%
    summarise(n = n(), n_pos = sum(seropos), .groups = "drop")
  ci_mat <- t(mapply(binom_ci, tab$n_pos, tab$n))
  
  all_obs[[ag]] <- data.frame(
    antigen = ag, age_mid = tab$age_mid, n = tab$n, n_pos = tab$n_pos,
    pct   = tab$n_pos / tab$n,
    lower = ci_mat[, "lower"], upper = ci_mat[, "upper"]
  )
  
  prob_sero <- ps$prob_seropositive
  pred_curve <- data.frame(
    agen  = stan_data$age,
    fit   = colMeans(prob_sero),
    lower = apply(prob_sero, 2, quantile, 0.025),
    upper = apply(prob_sero, 2, quantile, 0.975)
  )
  pred_curve <- pred_curve[order(pred_curve$agen), ]
  pred_curve$antigen <- ag
  all_pred[[ag]] <- pred_curve
  
  rm(ps, stan_data, prob_sero); gc()
}

obs_df  <- bind_rows(all_obs)
pred_df <- bind_rows(all_pred)

obs_df$antigen  <- factor(obs_df$antigen, levels = antigens)
pred_df$antigen <- factor(pred_df$antigen, levels = antigens)

write.csv(obs_df, "paper/tables/seroprevalence_observed_vs_predicted_finebin_age.csv", row.names = FALSE)


fig <- ggplot() +
  
  geom_ribbon(data = pred_df, aes(x = agen, ymin = lower, ymax = upper), fill = "#377eb8", alpha = 0.2) +
  geom_line(data = pred_df, aes(x = agen, y = fit, colour = "Predicted (model)"), linewidth = 0.9) +
  geom_errorbar(data = obs_df, aes(x = age_mid, ymin = lower, ymax = upper), width = 1.5, colour = "grey20") +
  geom_point(data = obs_df, aes(x = age_mid, y = pct, colour = "Observed (data-driven)"), size = 2) +
  facet_wrap(~ antigen, ncol = 4, labeller = as_labeller(ANTIGEN_LABELS)) +
  scale_colour_manual(name = NULL, values = c("Observed (data-driven)" = "grey20", "Predicted (model)" = "#377eb8")) +
  scale_y_continuous("Seroprevalence", labels = percent, limits = c(0, 1)) +
  labs(x = sprintf("Age (years)", BIN_WIDTH)) +
  theme_bw(base_size = 16) +
  theme(
    strip.text       = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "grey92", colour = NA),
    axis.title       = element_text(size = 16, face = "bold"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 13)
  )

ggsave("paper/seroprevalence_observed_vs_predicted_finebin_age.png", fig,
       width = 16, height = 10, dpi = 300, bg = "white")

