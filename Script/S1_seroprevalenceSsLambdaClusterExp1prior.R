# ============================================================================
# SEROPREVALENCE VS LAMBDA — CLUSTER LEVEL (4-comp rho-age, Exp(1) prior)
# ============================================================================
#
# Outputs:
#   paper/seroprevalence_vs_lambda_cluster_exp1prior_faceted.png
#   paper/tables/seroprevalence_vs_lambda_cluster_exp1prior_data.csv
#   paper/tables/seroprevalence_vs_lambda_cor_exp1prior.csv
# ============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(sn)   # dsn()



antigens     <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")
model_subdir <- "skew_normal_fixed_r8_1mu_rho_age_cluster_exp1prior"
out_dir      <- "paper"
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

mfi_cols <- setNames(paste0("log_mfi_", antigens), antigens)

# ============================================================================
# 1. MFI SEROPOS/SERONEG CROSSOVER THRESHOLD PER ANTIGEN
# ============================================================================

find_mfi_threshold <- function(mu_neg, sigma_neg, skew_neg,
                               mu_pos, sigma_pos, skew_pos) {
  f <- function(x) {
    dsn(x, xi = mu_neg, omega = sigma_neg, alpha = skew_neg) -
      dsn(x, xi = mu_pos, omega = sigma_pos, alpha = skew_pos)
  }
  tryCatch(
    uniroot(f, interval = c(mu_neg, mu_pos), tol = 1e-6)$root,
    error = function(e) (mu_neg + mu_pos) / 2
  )
}

threshold_table <- map_dfr(antigens, function(ag) {
  ps_file <- file.path("results", ag, model_subdir, "posterior_samples.rds")
  if (!file.exists(ps_file)) { warning("Missing posterior for ", ag); return(NULL) }
  ps  <- readRDS(ps_file)
  thr <- find_mfi_threshold(
    median(ps$mu_neg),   median(ps$sigma_neg), median(ps$skew_neg),
    median(ps$mu_pos),   median(ps$sigma_pos), median(ps$skew_pos)
  )
  data.frame(antigen = ag, threshold = thr)
})

print(threshold_table)

# ============================================================================
# 2. EMPIRICAL SEROPREVALENCE PER CLUSTER (same exclusion as model fitting)
# ============================================================================

data_clean <- read.csv("Data/data_cleaned.csv")

sero_long <- data_clean %>%
  filter(!is.na(cluster)) %>%
  filter(!(cluster %in% c(5, 32))) %>%
  select(cluster, all_of(unname(mfi_cols))) %>%
  pivot_longer(cols = all_of(unname(mfi_cols)),
               names_to = "mfi_col", values_to = "log_mfi") %>%
  mutate(antigen = names(mfi_cols)[match(mfi_col, mfi_cols)]) %>%
  filter(!is.na(log_mfi))

sero_prev <- threshold_table %>%
  left_join(sero_long, by = "antigen") %>%
  mutate(seropos = as.integer(log_mfi > threshold)) %>%
  group_by(antigen, cluster) %>%
  summarise(
    n                = n(),
    seroprev_pct     = mean(seropos) * 100,
    .groups = "drop"
  )

# ============================================================================
# 3. LAMBDA PER CLUSTER — Exp(1) PRIOR MODEL (from precomputed CSV)
# ============================================================================

load_lambda <- function(ag) {
  f <- file.path("results", ag, model_subdir, "lambda_per_cluster.csv")
  if (!file.exists(f)) { warning("Missing lambda_per_cluster.csv for ", ag); return(NULL) }
  lpc <- read.csv(f)
  data.frame(
    antigen       = ag,
    cluster       = as.integer(lpc$cluster),
    lambda_median = lpc$median,
    lambda_lo     = lpc$q025,
    lambda_hi     = lpc$q975
  )
}

lambda_all <- map_dfr(antigens, load_lambda)

# ============================================================================
# 4. JOIN
# ============================================================================

combined <- lambda_all %>%
  inner_join(sero_prev, by = c("antigen", "cluster")) %>%
  filter(!is.na(seroprev_pct))

write.csv(combined,
          file.path(out_dir, "tables/seroprevalence_vs_lambda_cluster_exp1prior_data.csv"),
          row.names = FALSE)

# ============================================================================
# 5. SCATTER PLOT — ALL ANTIGENS FACETED
# ============================================================================

cors <- combined %>%
  group_by(antigen) %>%
  summarise(
    r   = cor(seroprev_pct, lambda_median, method = "spearman", use = "complete.obs"),
    pv  = cor.test(seroprev_pct, lambda_median, method = "spearman")$p.value,
    x   = max(seroprev_pct, na.rm = TRUE),
    y   = max(lambda_hi,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("rho = %.2f\np = %.3f", r, pv))

write.csv(cors,
          file.path(out_dir, "tables/seroprevalence_vs_lambda_cor_exp1prior.csv"),
          row.names = FALSE)

p_facet <- ggplot(combined, aes(x = seroprev_pct, y = lambda_median)) +
  geom_errorbar(aes(ymin = lambda_lo, ymax = lambda_hi),
                width = 0, linewidth = 0.35, colour = "grey60", alpha = 0.6) +
  geom_point(colour = "#009E73", alpha = 0.85, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, colour = "#D55E00",
              fill = "#D55E00", alpha = 0.15, linewidth = 0.8) +
  geom_text(aes(label = cluster), vjust = -0.8, size = 2.4, colour = "grey30") +
  geom_text(data = cors, aes(x = x, y = y, label = label),
            hjust = 1, vjust = 1, size = 3.2, colour = "grey20", inherit.aes = FALSE) +
  facet_wrap(~ antigen, scales = "free_y", ncol = 4) +
  labs(
    title    = "Force of infection vs seroprevalence by village",
    subtitle = "4-compartment rho-age model, r fixed at 8/year, lambda ~ Exp(1) prior  |  Points = villages, error bars = 95% CrI",
    x        = "Seroprevalence (%)",
    y        = expression(paste("Seroconversion rate ", lambda, " (per year)"))
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, colour = "grey40"),
    strip.text       = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95", colour = NA),
    axis.title       = element_text(face = "bold", size = 11),
    legend.position  = "bottom"
  )

ggsave(file.path(out_dir, "seroprevalence_vs_lambda_cluster_exp1prior_faceted.png"),
       plot = p_facet, width = 16, height = 10, dpi = 300)


