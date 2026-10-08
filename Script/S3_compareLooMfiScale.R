################################################################################
# LOO COMPARISON — MFI-COMPARABLE SCALE
#
# Models compared (all on the same MFI-only log-likelihood scale):
#   M1: 2-compartment, constant rho       (2compartments_skewnormal)
#   M2: 2-compartment, age-varying rho    (skew_normal_2comp_rho_age)
#   M3: 4-compartment, age-varying rho    (skew_normal_fixed_r8_1mu_rho_age)
#   M4: 4-compartment, constant rho       (skew_normal_fixed_r8)
#
# The 2-comp models use log p(MFI | age, θ):
#   log[ (1 − seroprevalence_i) × f_neg(MFI_i) + seroprevalence_i × f_pos(MFI_i) ]
#
# The 4-comp model stores log p(MFI, PCR | age, θ) — a joint likelihood that
# conditions on the observed PCR status.
#
# Recompute the marginal MFI log-likelihood for the 4-comp
# model directly from stored posteriors:
#
#   ll_mfi[d, i] = log[ (1 − prob_seropositive[d,i]) × f_neg(MFI_i; θ_d)
#                       + prob_seropositive[d,i]      × f_pos(MFI_i; θ_d) ]
#
# This is identical in form to the 2-comp likelihood — both mix two skew-normal
# components weighted by the posterior seroprevalence — and is the proper
# MFI-only scale for LOO comparison.
#
# Output:
#   results/loo_comparison_mfi_scale.csv
#   results/loo_comparison_mfi_scale_wide.csv
#   results/parameter_table_all_models_raw.csv
#   results/parameter_table_all_models_wide.csv
#   results/loo_mfi_scale_elpd.png
#   results/loo_mfi_scale_delta.png
#   results/loo_mfi_scale_params.png
################################################################################

rm(list = ls())


suppressPackageStartupMessages({
  library(loo)
  library(sn)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

antigens <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")

model_dirs <- c(
  "2comp_const"   = "2compartments_skewnormal",
  "2comp_rho_age" = "skew_normal_2comp_rho_age",
  "4comp_const"   = "skew_normal_fixed_r8",
  "4comp_rho_age" = "skew_normal_fixed_r8_1mu_rho_age"
)

model_labels <- c(
  "2comp_const"   = "2-comp, constant rho",
  "2comp_rho_age" = "2-comp, age-varying rho",
  "4comp_const"   = "4-comp, constant rho",
  "4comp_rho_age" = "4-comp, age-varying rho"
)

################################################################################
# HELPER: summarise a vector of posterior draws
################################################################################

summ <- function(x) {
  if (is.null(x) || all(is.na(x))) return(c(mean=NA, sd=NA, q025=NA, q50=NA, q975=NA))
  q <- unname(quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  c(mean = mean(x, na.rm = TRUE),
    sd   = sd(x,   na.rm = TRUE),
    q025 = q[1],
    q50  = q[2],
    q975 = q[3])
}

fmt <- function(x, digits=3) {
  if (all(is.na(x))) return(NA_character_)
  s <- summ(x)
  sprintf("%.*f [%.*f, %.*f]",
          digits, s["mean"], digits, s["q025"], digits, s["q975"])
}

################################################################################
# EXTRACT PARAMETERS — harmonised across models
################################################################################

extract_params <- function(ps, model_key) {
  
  rows <- list()
  
  add <- function(par_name, draws) {
    if (is.null(draws) || all(is.na(draws))) return()
    s <- summ(draws)
    rows[[length(rows)+1]] <<- data.frame(
      parameter = par_name,
      mean  = s["mean"],  sd   = s["sd"],
      q025  = s["q025"],  q50  = s["q50"], q975 = s["q975"],
      stringsAsFactors = FALSE
    )
  }
  
  # ── lambda ──────────────────────────────────────────────────────────────────
  add("lambda", ps$lambda)
  
  # ── alpha (4-comp only) ─────────────────────────────────────────────────────
  add("alpha", ps$alpha)
  
  # ── rho ─────────────────────────────────────────────────────────────────────
  if (model_key %in% c("2comp_const", "4comp_const")) {
    add("rho", ps$rho)
  } else {
    for (k in 1:3) add(sprintf("rho[%d]", k), ps$rho[, k])
  }
  
  # ── MFI parameters — harmonise names across models ─────────────────────────
  # M1 (2comp_const): old naming mu_0/mu_1/epsilon_0/epsilon_1/alpha_0/alpha_1
  # M4 (4comp_const): uses delta_mu (not delta_sero) for mu_pos - mu_neg offset;
  #                   delta_mu is not extracted in the run script so ps$delta_mu
  #                   is NULL in saved RDS — mu_pos is reported directly instead
  if (model_key == "2comp_const") {
    add("mu_neg",   ps$mu_0)
    add("sigma_neg",ps$epsilon_0)
    add("skew_neg", ps$alpha_0)
    add("mu_pos",   ps$mu_0 + ps$mu_1)
    add("sigma_pos",ps$epsilon_1)
    add("skew_pos", ps$alpha_1)
  } else if (model_key == "4comp_const") {
    add("mu_neg",    ps$mu_neg)
    add("sigma_neg", ps$sigma_neg)
    add("skew_neg",  ps$skew_neg)
    add("delta_mu",  ps$delta_mu)   # NULL in existing RDS; mu_pos shown directly
    add("mu_pos",    ps$mu_pos)
    add("sigma_pos", ps$sigma_pos)
    add("skew_pos",  ps$skew_pos)
  } else {
    add("mu_neg",    ps$mu_neg)
    add("sigma_neg", ps$sigma_neg)
    add("skew_neg",  ps$skew_neg)
    add("delta_sero",ps$delta_sero)
    add("mu_pos",    ps$mu_pos)
    add("sigma_pos", ps$sigma_pos)
    add("skew_pos",  ps$skew_pos)
  }
  
  do.call(rbind, rows)
}

################################################################################
# MAIN LOOP
################################################################################

loo_rows   <- list()
param_rows <- list()

for (ag in antigens) {
  
  # Load stan data for MFI vector — shared across all 4-comp models
  sd4_path <- sprintf("results/%s/%s/stan_data.rds", ag, model_dirs["4comp_rho_age"])
  if (!file.exists(sd4_path))
    sd4_path <- sprintf("results/%s/%s/stan_data.rds", ag, model_dirs["4comp_const"])
  sd4 <- readRDS(sd4_path)
  log_mfi_vec <- sd4$log_mfi
  
  loo_list <- list()
  
  for (mkey in names(model_dirs)) {
    
    path_ps <- sprintf("results/%s/%s/posterior_samples.rds",
                       ag, model_dirs[mkey])
    
    if (!file.exists(path_ps)) {
      next
    }
    
    ps <- readRDS(path_ps)
    
    # ── log-likelihood matrix (draws × N) ─────────────────────────────────────
    ll_mat <- ps$log_lik  # already [draws × N]
    
    # ── 4-comp: recompute marginal MFI log-likelihood ────────────────────────
    # Subtracting log p(PCR | θ) gives the CONDITIONAL likelihood
    # log p(MFI | PCR_observed, θ), which is biased upward vs the 2-comp models.
    # Instead, recompute log p(MFI | age, θ) directly as a skew-normal mixture
    # weighted by prob_seropositive — identical in form to the 2-comp likelihood.
    if (mkey %in% c("4comp_const", "4comp_rho_age")) {
      seroprev_mat <- ps$prob_seropositive  # [draws × N]
      n_d <- nrow(seroprev_mat)
      ll_mat <- do.call(rbind, lapply(seq_len(n_d), function(d) {
        f_neg <- dsn(log_mfi_vec,
                     xi    = ps$mu_neg[d],
                     omega = ps$sigma_neg[d],
                     alpha = ps$skew_neg[d])
        f_pos <- dsn(log_mfi_vec,
                     xi    = ps$mu_pos[d],
                     omega = ps$sigma_pos[d],
                     alpha = ps$skew_pos[d])
        sp <- pmin(pmax(seroprev_mat[d, ], 0), 1)
        log(pmax((1 - sp) * f_neg + sp * f_pos, 1e-300))
      }))
    }
    
    # ── LOO ──────────────────────────────────────────────────────────────────
    loo_res <- tryCatch(
      loo(ll_mat, cores = 4),
      error = function(e) NULL
    )
    
    if (!is.null(loo_res)) {
      elpd  <- loo_res$estimates["elpd_loo", "Estimate"]
      elpd_se <- loo_res$estimates["elpd_loo", "SE"]
      p_loo <- loo_res$estimates["p_loo", "Estimate"]
      looic <- loo_res$estimates["looic", "Estimate"]
      loo_list[[mkey]] <- loo_res
      
      loo_rows[[length(loo_rows)+1]] <- data.frame(
        antigen   = ag,
        model_key = mkey,
        model     = model_labels[mkey],
        elpd      = elpd,
        elpd_se   = elpd_se,
        p_loo     = p_loo,
        looic     = looic,
        stringsAsFactors = FALSE
      )
    }
    
    # ── Parameters ──────────────────────────────────────────────────────────
    pars <- extract_params(ps, mkey)
    pars$antigen   <- ag
    pars$model_key <- mkey
    pars$model     <- model_labels[mkey]
    param_rows[[length(param_rows)+1]] <- pars
  }
  
  # ── Delta ELPD within antigen ────────────────────────────────────────────
  if (length(loo_list) >= 2) {
    loo_names <- names(loo_list)
    best_elpd <- max(sapply(loo_list, function(l)
      l$estimates["elpd_loo", "Estimate"]))
    for (mkey in loo_names) {
      elpd_i <- loo_list[[mkey]]$estimates["elpd_loo", "Estimate"]
      delta  <- elpd_i - best_elpd
      idx <- which(sapply(loo_rows, function(r)
        r$antigen == ag & r$model_key == mkey))
      if (length(idx) == 1)
        loo_rows[[idx]]$delta_elpd <- delta
    }
  }
}

################################################################################
# ASSEMBLE TABLES
################################################################################

loo_df   <- bind_rows(loo_rows)
param_df <- bind_rows(param_rows)

# Fill missing delta_elpd (best model = 0)
if (!"delta_elpd" %in% names(loo_df)) loo_df$delta_elpd <- NA_real_
loo_df$delta_elpd[is.na(loo_df$delta_elpd)] <- 0

# ── LOO wide table (antigens × models) ──────────────────────────────────────
loo_wide <- loo_df %>%
  mutate(cell = sprintf("%.1f (±%.1f)", elpd, elpd_se)) %>%
  select(antigen, model, cell) %>%
  pivot_wider(names_from = model, values_from = cell)

print(loo_wide, n=Inf)

# ── Delta ELPD table ─────────────────────────────────────────────────────────
delta_wide <- loo_df %>%
  mutate(cell = sprintf("%.1f", delta_elpd)) %>%
  select(antigen, model, cell) %>%
  pivot_wider(names_from = model, values_from = cell)

print(delta_wide, n=Inf)

################################################################################
# PARAMETER TABLE — formatted wide
################################################################################

param_order <- c("lambda", "alpha",
                 "rho", "rho[1]", "rho[2]", "rho[3]",
                 "mu_neg", "sigma_neg", "skew_neg",
                 "delta_mu", "delta_sero",
                 "mu_pos", "sigma_pos", "skew_pos")

param_fmt <- param_df %>%
  mutate(
    value_fmt = sprintf("%.3f [%.3f, %.3f]", mean, q025, q975),
    parameter = factor(parameter, levels = param_order)
  ) %>%
  arrange(antigen, model_key, parameter) %>%
  select(antigen, model, parameter, value_fmt) %>%
  pivot_wider(names_from = model, values_from = value_fmt)

print(param_fmt, n=Inf)

################################################################################
# SAVE CSV
################################################################################

write.csv(loo_df,    "results/loo_comparison_mfi_scale.csv",    row.names=FALSE)
write.csv(loo_wide,  "results/loo_comparison_mfi_scale_wide.csv", row.names=FALSE)
write.csv(param_df,  "results/parameter_table_all_models_raw.csv",  row.names=FALSE)
write.csv(param_fmt, "results/parameter_table_all_models_wide.csv", row.names=FALSE)

################################################################################
# PLOT 1 — ELPD by model and antigen
################################################################################

model_order <- c("2-comp, constant rho", "2-comp, age-varying rho",
                 "4-comp, constant rho", "4-comp, age-varying rho")
model_colors <- c("2-comp, constant rho"   = "#1b9e77",
                  "2-comp, age-varying rho" = "#7570b3",
                  "4-comp, constant rho"    = "#e7298a",
                  "4-comp, age-varying rho" = "#d95f02")

p_elpd <- loo_df %>%
  mutate(model = factor(model, levels = model_order)) %>%
  ggplot(aes(x = antigen, y = elpd, colour = model, group = model)) +
  geom_errorbar(aes(ymin = elpd - elpd_se, ymax = elpd + elpd_se),
                width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_line(position = position_dodge(0.5), linewidth = 0.5, linetype = "dashed") +
  scale_colour_manual(values = model_colors, name = "Model") +
  labs(
    title   = "LOO-ELPD comparison — MFI-only scale",
    subtitle= "4-comp marginal MFI log-lik recomputed from prob_seropositive × skew-normal mixture",
    x = "Antigen", y = "ELPD (higher = better)"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "grey40"))

ggsave("results/loo_mfi_scale_elpd.png", p_elpd, width = 11, height = 6, dpi = 300)
ggsave("results/loo_mfi_scale_elpd.pdf", p_elpd, width = 11, height = 6)

################################################################################
# PLOT 2 — Delta ELPD (relative to best model per antigen)
################################################################################

p_delta <- loo_df %>%
  mutate(model = factor(model, levels = model_order)) %>%
  ggplot(aes(x = antigen, y = delta_elpd, colour = model, group = model)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = c(-4, -10), colour = "grey70",
             linewidth = 0.4, linetype = "dotted") +
  geom_errorbar(aes(ymin = delta_elpd - elpd_se, ymax = delta_elpd + elpd_se),
                width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_colour_manual(values = model_colors, name = "Model") +
  scale_y_continuous("delta-ELPD vs best model (0 = best)",
                     breaks = pretty(loo_df$delta_elpd, n=8)) +
  labs(
    title    = "delta-ELPD — MFI-only scale (comparable across 2- and 4-comp)",
    subtitle = "Dotted lines: −4 and −10 ELPD thresholds",
    x = "Antigen"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "grey40"))

ggsave("results/loo_mfi_scale_delta.png", p_delta, width = 11, height = 6, dpi = 300)
ggsave("results/loo_mfi_scale_delta.pdf", p_delta, width = 11, height = 6)

################################################################################
# PLOT 3 — Key parameters across models and antigens
################################################################################

key_params <- c("lambda", "rho[1]", "rho[2]", "rho[3]", "rho")

p_params <- param_df %>%
  filter(parameter %in% key_params) %>%
  mutate(
    parameter = recode(parameter, "rho" = "rho (constant)"),
    model = factor(model, levels = model_order)
  ) %>%
  ggplot(aes(x = antigen, y = mean, colour = model,
             ymin = q025, ymax = q975, group = model)) +
  geom_errorbar(width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = model_colors, name = "Model") +
  labs(title = "Parameter estimates across models and antigens",
       subtitle = "Mean ± 95% CrI",
       x = "Antigen", y = "Posterior mean [95% CrI]") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

ggsave("results/loo_mfi_scale_params.png", p_params,
       width = 12, height = 9, dpi = 300)
ggsave("results/loo_mfi_scale_params.pdf", p_params, width = 12, height = 9)
