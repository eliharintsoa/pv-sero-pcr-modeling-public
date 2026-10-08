# =============================================================================
# SUPPLEMENTARY — Village (cluster)-level lambda caterpillar plot, all 8 antigens
#   4-compartment, age-varying rho, Exp(1)-prior cluster model (production model)
#
# y-axis = clusters (31 villages) ordered by median lambda (consistent
# ordering derived from mean-of-medians across antigens). Point = posterior
# median; bar = 95% credible interval.
#
# Output: paper/figure6_caterpillar_lambda_village.png / .pdf
# =============================================================================



suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

dir.create("paper", showWarnings = FALSE)

ANTIGENS <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")

# ---------------------------------------------------------------------------
# Build tidy data frame: one row per cluster × antigen
# ---------------------------------------------------------------------------
extract_cluster_lambda <- function(model_dir_template, fallback_ids = NULL) {
  do.call(rbind, lapply(ANTIGENS, function(ag) {
    dir  <- sprintf(model_dir_template, ag)
    ps   <- readRDS(file.path(dir, "posterior_samples.rds"))
    sd   <- readRDS(file.path(dir, "stan_data.rds"))
    
    draws <- ps$lambda                    # matrix [n_draws x N_groups]
    n_cl  <- ncol(draws)
    names <- if (!is.null(sd$cluster_names) && length(sd$cluster_names) == n_cl)
      paste0("Cl.", sd$cluster_names)
    else if (!is.null(fallback_ids) && length(fallback_ids) == n_cl)
      paste0("Cl.", fallback_ids)
    else paste0("Cl.", seq_len(n_cl))
    
    data.frame(
      Antigen = ag,
      Cluster = names,
      median  = apply(draws, 2, median),
      lo95    = apply(draws, 2, quantile, 0.025),
      hi95    = apply(draws, 2, quantile, 0.975),
      lo50    = apply(draws, 2, quantile, 0.25),
      hi50    = apply(draws, 2, quantile, 0.75),
      stringsAsFactors = FALSE
    )
  }))
}

# ---------------------------------------------------------------------------
# Shared cluster ordering (mean-of-medians across antigens)
# so the same cluster sits at the same rank across facets
# ---------------------------------------------------------------------------
cluster_order <- function(df) {
  df |>
    group_by(Cluster) |>
    summarise(mean_med = mean(median), .groups = "drop") |>
    arrange(mean_med) |>
    pull(Cluster)
}

# ---------------------------------------------------------------------------
# Plotting function
# ---------------------------------------------------------------------------
caterpillar_cluster_plot <- function(df, title, colour) {
  ord <- cluster_order(df)
  df  <- df |>
    mutate(
      Antigen = factor(Antigen, levels = ANTIGENS),
      Cluster = factor(Cluster, levels = ord)
    )
  
  ggplot(df, aes(y = Cluster, x = median)) +
    # 95% CrI
    geom_segment(aes(x = lo95, xend = hi95, yend = Cluster),
                 colour = colour, linewidth = 0.45, alpha = 0.45) +
    # median point
    geom_point(colour = colour, size = 1.4) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    facet_wrap(~ Antigen, nrow = 2, scales = "free_x") +
    labs(
      title = title,
      x     = expression(lambda ~ "(seroconversion rate, year"^{-1}*")"),
      y     = "Cluster"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
      strip.text      = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "grey95", colour = NA),
      axis.text.y     = element_text(size = 7),
      axis.text.x     = element_text(size = 7),
      axis.title      = element_text(size = 9),
      panel.spacing.x = unit(0.6, "lines"),
      panel.spacing.y = unit(0.4, "lines")
    )
}

# ---------------------------------------------------------------------------
# 4-compartment, age-varying rho, Exp(1)-prior cluster model
# ---------------------------------------------------------------------------
df_village <- extract_cluster_lambda(
  "results/%s/skew_normal_fixed_r8_1mu_rho_age_cluster_exp1prior",
  fallback_ids = setdiff(1:33, c(5, 32))
)

p_village <- caterpillar_cluster_plot(
  df_village,
  title  = "Village-level seroconversion rate (lambda), 4-compartment rho-age model",
  colour = "#0072B2"
)

ggsave("paper/figure6_caterpillar_lambda_village.png", p_village,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave("paper/figure6_caterpillar_lambda_village.pdf", p_village,
       width = 14, height = 10)
