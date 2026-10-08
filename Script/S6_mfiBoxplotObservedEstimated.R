################################################################################
# MFI BOXPLOT — OBSERVED vs ESTIMATED, BY AGE
#
# Compares observed log-MFI (by 5-year age bin, stratified by PCR status)
# against posterior-predictive draws from four model variants:
#   - 4c_fixed_rho : 4-comp model, constant rho
#   - 4c_rho_age   : 4-comp model, age-varying rho (1mu)
#   - 2c_fixed_rho : 2-comp model, constant rho
#   - 2c_rho_age   : 2-comp model, age-varying rho
#
# For each model, 100 posterior draws are used to simulate individual-level
# log-MFI values (same ages/PCR status as the observed data), giving a
# posterior-predictive distribution directly comparable to the observed one.
# 2-comp models have no PCR compartment; their predicted distributions use
# the age-specific marginal seroprevalence mixture without PCR conditioning.
#
# Output: results/{antigen}/plots/mfi_boxplot_observed_estimated.png/.pdf
################################################################################


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(sn)
})

# compartment_model_2comp.R sources compartment_model.R internally, so both
# 4-comp and 2-comp functions are available after this single source call.
source("compartment_model_2comp.R")

r_fixed     <- 8
dir_4c_fr   <- sprintf("skew_normal_fixed_r%.0f",             r_fixed)
dir_4c_ra   <- sprintf("skew_normal_fixed_r%.0f_1mu_rho_age", r_fixed)
dir_2c_fr   <- "2compartments_skewnormal"
dir_2c_ra   <- "skew_normal_2comp_rho_age"
n_draws_box <- 100

antigens <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")
ANTIGEN_LABELS <- c(EBP = "EBP", LF005 = "Pv-fam-a", LF010 = "MSP5", LF016 = "MSP1-19",
                    MSP8 = "MSP8", PTEX150 = "PTEX150", PvCSS = "PvCSS",
                    RBP2b.P87 = "RBP2b")

col_obs   <- "grey70"
col_4c_fr <- "#2166ac"
col_4c_ra <- "#d6604d"
col_2c_fr <- "lightblue"
col_2c_ra <- "#984ea3"

model_keys <- c("4c_fixed_rho", "4c_rho_age", "2c_fixed_rho", "2c_rho_age")

model_colours <- c(
  "4c_fixed_rho" = col_4c_fr,
  "4c_rho_age"   = col_4c_ra,
  "2c_fixed_rho" = col_2c_fr,
  "2c_rho_age"   = col_2c_ra
)

model_labels_map <- c(
  "4c_fixed_rho" = sprintf("4-comp const-rho (r=%g)", r_fixed),
  "4c_rho_age"   = sprintf("4-comp rho-age (r=%g)",   r_fixed),
  "2c_fixed_rho" = "2-comp const-rho",
  "2c_rho_age"   = "2-comp rho-age"
)

# Normalise MFI parameter names: 2-comp const-rho uses old naming convention
# (mu_0/mu_1/epsilon_0/epsilon_1/alpha_0/alpha_1); all others use standard names.
norm_mfi_params <- function(ps) {
  list(
    mu_neg    = if (!is.null(ps$mu_neg))    ps$mu_neg    else ps$mu_0,
    sigma_neg = if (!is.null(ps$sigma_neg)) ps$sigma_neg else ps$epsilon_0,
    skew_neg  = if (!is.null(ps$skew_neg))  ps$skew_neg  else ps$alpha_0,
    mu_pos    = if (!is.null(ps$mu_pos))    ps$mu_pos    else ps$mu_0 + ps$mu_1,
    sigma_pos = if (!is.null(ps$sigma_pos)) ps$sigma_pos else ps$epsilon_1,
    skew_pos  = if (!is.null(ps$skew_pos))  ps$skew_pos  else ps$alpha_1
  )
}

# Prefer best-chains file when available, fall back to full posterior.
load_posterior <- function(ag, subdir) {
  best <- file.path("results", ag, subdir, "posterior_samples_best_chains.rds")
  full <- file.path("results", ag, subdir, "posterior_samples.rds")
  if (file.exists(best)) best else if (file.exists(full)) full else NULL
}

for (ag in antigens) {
  
  data_file <- sprintf("Data/stan_data_%s.rds", ag)
  if (!file.exists(data_file)) next
  sd <- readRDS(data_file)
  
  ages    <- sd$age
  log_mfi <- sd$log_mfi
  pcr     <- sd$pcr
  n       <- sd$N
  
  n_neg <- sum(pcr == 0)
  n_pos <- sum(pcr == 1)
  
  age_breaks <- seq(0, max(ages) + 5, by = 5)
  age_groups <- cut(ages, breaks = age_breaks, include.lowest = TRUE,
                    labels = head(age_breaks, -1) + diff(age_breaks) / 2)
  
  out_dir <- file.path("results", ag, "plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── Model 1: 4-comp, constant rho ───────────────────────────────────────────
  path_4c_fr <- load_posterior(ag, dir_4c_fr)
  has_4c_fr  <- !is.null(path_4c_fr)
  if (has_4c_fr) {
    ps_m1 <- readRDS(path_4c_fr)
    p_m1  <- norm_mfi_params(ps_m1)
    set.seed(123)
    draw_idx  <- sample(length(ps_m1$lambda), min(n_draws_box, length(ps_m1$lambda)))
    df_sim_m1 <- do.call(rbind, lapply(draw_idx, function(d) {
      s_d <- sapply(ages, function(a)
        calculate_compartments(a, ps_m1$lambda[d], ps_m1$alpha[d], ps_m1$rho[d], r_fixed))
      data.frame(
        log_mfi   = simulate_mfi(s_d, pcr,
                                 p_m1$mu_neg[d],   p_m1$sigma_neg[d], p_m1$skew_neg[d],
                                 p_m1$mu_pos[d],   p_m1$sigma_pos[d], p_m1$skew_pos[d]),
        age_bin   = age_groups,
        pcr_label = ifelse(pcr == 0, "PCR negative", "PCR positive"),
        Source    = "4c_fixed_rho"
      )
    }))
  }
  
  # ── Model 2: 4-comp, rho-age ────────────────────────────────────────────────
  path_4c_ra <- load_posterior(ag, dir_4c_ra)
  has_4c_ra  <- !is.null(path_4c_ra)
  if (has_4c_ra) {
    ps_m2 <- readRDS(path_4c_ra)
    if (is.null(dim(ps_m2$rho)))
      ps_m2$rho <- cbind(ps_m2[["rho[1]"]], ps_m2[["rho[2]"]], ps_m2[["rho[3]"]])
    p_m2 <- norm_mfi_params(ps_m2)
    set.seed(123)
    draw_idx  <- sample(length(ps_m2$lambda), min(n_draws_box, length(ps_m2$lambda)))
    df_sim_m2 <- do.call(rbind, lapply(draw_idx, function(d) {
      rho_d <- ps_m2$rho[d, ]
      s_d   <- sapply(ages, function(a)
        calculate_compartments_age_varying(a, ps_m2$lambda[d], ps_m2$alpha[d], rho_d, r_fixed))
      data.frame(
        log_mfi   = simulate_mfi(s_d, pcr,
                                 p_m2$mu_neg[d], p_m2$sigma_neg[d], p_m2$skew_neg[d],
                                 p_m2$mu_pos[d], p_m2$sigma_pos[d], p_m2$skew_pos[d]),
        age_bin   = age_groups,
        pcr_label = ifelse(pcr == 0, "PCR negative", "PCR positive"),
        Source    = "4c_rho_age"
      )
    }))
  }
  
  # ── Model 3: 2-comp, constant rho ───────────────────────────────────────────
  path_2c_fr <- load_posterior(ag, dir_2c_fr)
  has_2c_fr  <- !is.null(path_2c_fr)
  if (has_2c_fr) {
    ps_m3 <- readRDS(path_2c_fr)
    p_m3  <- norm_mfi_params(ps_m3)
    set.seed(123)
    draw_idx  <- sample(length(ps_m3$lambda), min(n_draws_box, length(ps_m3$lambda)))
    df_sim_m3 <- do.call(rbind, lapply(draw_idx, function(d) {
      s_d <- sapply(ages, function(a)
        calculate_2comp(a, ps_m3$lambda[d], ps_m3$rho[d]))
      data.frame(
        log_mfi   = simulate_2comp_mfi(s_d,
                                       p_m3$mu_neg[d], p_m3$sigma_neg[d], p_m3$skew_neg[d],
                                       p_m3$mu_pos[d], p_m3$sigma_pos[d], p_m3$skew_pos[d]),
        age_bin   = age_groups,
        pcr_label = ifelse(pcr == 0, "PCR negative", "PCR positive"),
        Source    = "2c_fixed_rho"
      )
    }))
  }
  
  # ── Model 4: 2-comp, rho-age ────────────────────────────────────────────────
  path_2c_ra <- load_posterior(ag, dir_2c_ra)
  has_2c_ra  <- !is.null(path_2c_ra)
  if (has_2c_ra) {
    ps_m4 <- readRDS(path_2c_ra)
    if (is.null(dim(ps_m4$rho)))
      ps_m4$rho <- cbind(ps_m4[["rho[1]"]], ps_m4[["rho[2]"]], ps_m4[["rho[3]"]])
    p_m4 <- norm_mfi_params(ps_m4)
    set.seed(123)
    draw_idx  <- sample(length(ps_m4$lambda), min(n_draws_box, length(ps_m4$lambda)))
    df_sim_m4 <- do.call(rbind, lapply(draw_idx, function(d) {
      rho_d <- ps_m4$rho[d, ]
      s_d   <- sapply(ages, function(a)
        calculate_compartments_age_varying_2comp(a, ps_m4$lambda[d], rho_d))
      data.frame(
        log_mfi   = simulate_2comp_mfi(s_d,
                                       p_m4$mu_neg[d], p_m4$sigma_neg[d], p_m4$skew_neg[d],
                                       p_m4$mu_pos[d], p_m4$sigma_pos[d], p_m4$skew_pos[d]),
        age_bin   = age_groups,
        pcr_label = ifelse(pcr == 0, "PCR negative", "PCR positive"),
        Source    = "2c_rho_age"
      )
    }))
  }
  
  if (!has_4c_fr && !has_4c_ra && !has_2c_fr && !has_2c_ra) next
  
  active_keys <- intersect(model_keys,
                           c(if (has_4c_fr) "4c_fixed_rho",
                             if (has_4c_ra) "4c_rho_age",
                             if (has_2c_fr) "2c_fixed_rho",
                             if (has_2c_ra) "2c_rho_age"))
  
  # ── Boxplot by 5-yr age bin, stratified by PCR status ──────────────────────
  df_obs_box <- data.frame(
    log_mfi   = log_mfi,
    age_bin   = age_groups,
    pcr_label = ifelse(pcr == 0, "PCR negative", "PCR positive"),
    Source    = "Observed"
  )
  df_box_all <- df_obs_box
  if (has_4c_fr) df_box_all <- bind_rows(df_box_all, df_sim_m1[!is.na(df_sim_m1$age_bin), ])
  if (has_4c_ra) df_box_all <- bind_rows(df_box_all, df_sim_m2[!is.na(df_sim_m2$age_bin), ])
  if (has_2c_fr) df_box_all <- bind_rows(df_box_all, df_sim_m3[!is.na(df_sim_m3$age_bin), ])
  if (has_2c_ra) df_box_all <- bind_rows(df_box_all, df_sim_m4[!is.na(df_sim_m4$age_bin), ])
  df_box_all <- df_box_all[!is.na(df_box_all$age_bin), ]
  
  source_levels  <- c("Observed", active_keys)
  source_colours <- c("Observed" = col_obs, model_colours)
  source_labels  <- c("Observed" = "Observed", model_labels_map)
  
  df_box_all$Source    <- factor(df_box_all$Source,    levels = source_levels)
  df_box_all$pcr_label <- factor(df_box_all$pcr_label,
                                 levels = c("PCR negative", "PCR positive"))
  
  p_box <- ggplot(df_box_all,
                  aes(x = age_bin, y = log_mfi, fill = Source, colour = Source)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4,
                 alpha = 0.55, linewidth = 0.35) +
    scale_fill_manual(values   = source_colours, labels = source_labels) +
    scale_colour_manual(values = source_colours, labels = source_labels) +
    facet_wrap(~ pcr_label, ncol = 1, scales = "free_y",
               labeller = labeller(pcr_label = c(
                 "PCR negative" = sprintf("PCR negative  (n = %d)", n_neg),
                 "PCR positive" = sprintf("PCR positive  (n = %d)", n_pos)
               ))) +
    labs(
      title  = sprintf("MFI boxplot — observed vs estimated — %s  (r = %g fixed)", ANTIGEN_LABELS[[ag]], r_fixed),
      x      = "Age (years)",
      y      = sprintf("log(MFI) — %s", ANTIGEN_LABELS[[ag]]),
      fill   = NULL,
      colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      strip.text      = element_text(face = "bold", size = 11)
    )
  
  ggsave(file.path(out_dir, "mfi_boxplot_observed_estimated.png"),
         p_box, width = 14, height = 10, dpi = 300)
  ggsave(file.path(out_dir, "mfi_boxplot_observed_estimated.pdf"),
         p_box, width = 14, height = 10)
}
