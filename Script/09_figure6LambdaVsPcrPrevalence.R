# ============================================================================
# FIGURE 6 — Madagascar: village-level force of infection (lambda) vs
# observed PCR prevalence, 4-compartment rho-age model, Exp(1)-prior cluster
# model
#
# Outputs:
#   paper/tables/spearman_pcr_vs_lambda_exp1prior.csv
#   paper/pcr_vs_lambda_cluster_4c_rho_age_exp1prior_faceted.png   (Figure 5)
#   paper/pcr_vs_lambda_cluster_4comp_rho_age_exp1prior_data.csv
#   paper/tables/lambda_range_by_antigen_4c_rhoage_exp1prior.csv
# ============================================================================

library(dplyr)
library(ggplot2)
library(purrr)

# Run from the repository root (Pv_sero_pcr_modeling/)

antigens <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")
ANTIGEN_LABELS <- c(
  "EBP"       = "EBP",
  "LF005"     = "Pv-fam-a",
  "LF010"     = "MSP5",
  "LF016"     = "MSP1-19",
  "MSP8"      = "MSP8",
  "PTEX150"   = "PTEX150",
  "PvCSS"     = "PvCSS",
  "RBP2b.P87" = "RBP2b"
)
out_dir  <- "paper"
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. PCR PREVALENCE BY CLUSTER
# ============================================================================
data_clean <- read.csv("Data/data_cleaned.csv")

pcr_prev <- data_clean %>%
  filter(!is.na(pcr), !is.na(cluster)) %>%
  filter(!(cluster %in% c(5, 32))) %>%
  mutate(pcr_pos = as.integer(pcr == "Pv infection")) %>%
  group_by(cluster) %>%
  summarise(n = n(), pcr_prev_pct = mean(pcr_pos, na.rm = TRUE) * 100, .groups = "drop")


# ============================================================================
# 2. LOAD LAMBDA — 4C RHO-AGE, EXP(1)-PRIOR CLUSTER MODEL
# ============================================================================
load_lambda_model <- function(ag, subdir) {
  ps_file <- file.path("results", ag, subdir, "posterior_samples.rds")
  sd_file <- file.path("results", ag, subdir, "stan_data.rds")
  if (!file.exists(ps_file) || !file.exists(sd_file)) {
    warning("Missing files for ", ag, " / ", subdir); return(NULL)
  }
  ps       <- readRDS(ps_file)
  sd       <- readRDS(sd_file)
  lam_mat  <- ps$lambda
  orig_ids <- setdiff(1:33, c(5, 32))
  data.frame(
    antigen       = ag,
    cluster       = as.integer(orig_ids),
    lambda_median = apply(lam_mat, 2, median),
    lambda_lo     = apply(lam_mat, 2, quantile, 0.025),
    lambda_med    = apply(lam_mat, 2, quantile, 0.500),
    lambda_hi     = apply(lam_mat, 2, quantile, 0.975),
    row.names     = NULL
  )
}

lambda_all_4c_ra_exp1 <- map_dfr(antigens, load_lambda_model,
                                 subdir = "skew_normal_fixed_r8_1mu_rho_age_cluster_exp1prior")


combined_4c_ra_exp1 <- lambda_all_4c_ra_exp1 %>%
  left_join(pcr_prev, by = "cluster") %>%
  filter(!is.na(pcr_prev_pct))

# ============================================================================
# 3. SPEARMAN CORRELATION PER ANTIGEN
# ============================================================================
cors_4c_ra_exp1 <- combined_4c_ra_exp1 %>%
  group_by(antigen) %>%
  summarise(
    r  = cor(pcr_prev_pct, lambda_med, method = "spearman", use = "complete.obs"),
    pv = cor.test(pcr_prev_pct, lambda_med, method = "spearman")$p.value,
    x  = max(pcr_prev_pct, na.rm = TRUE),
    y  = max(lambda_hi,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("rho = %.2f\np = %.3f", r, pv))

write.csv(cors_4c_ra_exp1[, c("antigen", "r", "pv")],
          file.path(out_dir, "tables", "spearman_pcr_vs_lambda_exp1prior.csv"),
          row.names = FALSE)

print(cors_4c_ra_exp1[, c("antigen", "r", "pv")])


# ============================================================================
# 4. FIGURE 7 — 4C rho-age, Exp(1) prior, faceted by antigen
# ============================================================================
p_4c_ra_exp1 <- ggplot(combined_4c_ra_exp1, aes(x = pcr_prev_pct, y = lambda_med)) +
  geom_errorbar(aes(ymin = lambda_lo, ymax = lambda_hi, linetype = "95% CrI"),
                width = 0, linewidth = 0.35, colour = "#0072B2", alpha = 0.5) +
  geom_point(aes(shape = "Cluster"), colour = "#0072B2", alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, colour = "#0072B2",
              fill = "#0072B2", alpha = 0.12, linewidth = 0.8) +
  geom_text(data = cors_4c_ra_exp1, aes(x = x, y = y, label = label),
            hjust = 1, vjust = 1, size = 5.2, colour = "grey20", inherit.aes = FALSE) +
  facet_wrap(~ antigen, scales = "free_y", ncol = 4, labeller = as_labeller(ANTIGEN_LABELS)) +
  scale_shape_manual(name = NULL, values = c("Cluster" = 16)) +
  scale_linetype_manual(name = NULL, values = c("95% CrI" = "solid")) +
  labs(
    title    = "Force of infection vs PCR prevalence by cluster",
    x        = "PCR prevalence (%)",
    y        = expression(paste("Seroconversion rate ", lambda, " (per year)"))
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    strip.text       = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey95", colour = NA),
    axis.title       = element_text(face = "bold", size = 14),
    legend.position  = "bottom"
  )

ggsave(file.path(out_dir, "pcr_vs_lambda_cluster_4c_rho_age_exp1prior_faceted.png"),
       plot = p_4c_ra_exp1, width = 16, height = 10, dpi = 300)

write.csv(combined_4c_ra_exp1,
          file.path(out_dir, "pcr_vs_lambda_cluster_4comp_rho_age_exp1prior_data.csv"),
          row.names = FALSE)

# ============================================================================
# 5. VILLAGE-LEVEL LAMBDA RANGE TABLE — EXP(1) PRIOR
# ============================================================================
lambda_range_exp1 <- combined_4c_ra_exp1 %>%
  group_by(antigen) %>%
  summarise(
    n_clusters        = n(),
    cluster_min       = cluster[which.min(lambda_med)],
    lambda_min_median = lambda_med[which.min(lambda_med)],
    lambda_min_lo     = lambda_lo[which.min(lambda_med)],
    lambda_min_hi     = lambda_hi[which.min(lambda_med)],
    cluster_max       = cluster[which.max(lambda_med)],
    lambda_max_median = lambda_med[which.max(lambda_med)],
    lambda_max_lo     = lambda_lo[which.max(lambda_med)],
    lambda_max_hi     = lambda_hi[which.max(lambda_med)],
    fold              = round(max(lambda_med) / min(lambda_med), 1),
    .groups = "drop"
  )

write.csv(lambda_range_exp1,
          file.path(out_dir, "tables", "lambda_range_by_antigen_4c_rhoage_exp1prior.csv"),
          row.names = FALSE)
