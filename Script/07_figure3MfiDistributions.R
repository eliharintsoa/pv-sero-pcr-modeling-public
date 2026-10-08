# =============================================================================
# FIGURE 4 — MFI distributions, observed vs fitted skew-normal mixture,
#            PCR-weighted (4C-rho-age model, all 8 antigens)
#
# Bars  = observed data (raw log-MFI counts, split by observed PCR status).
# Lines = model-predicted density: the fitted skew-normal mixture, weighted
#         by each PCR group's model-implied seropositivity proportion
#         (w_pos_pcrneg / w_pos_pcrpos). The total line is the PCR- line plus
#         the PCR+ contribution added on top, so its height is the total
#         across both PCR groups (matching how the bars are stacked).
#
# Output: paper/figure5_1mu_rho_age_hist_pcr_weighted.png  ( Figure 4)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(sn)
})

# Run from the repository root (Pv_sero_pcr_modeling/)
source("compartment_model.R")   # calculate_compartments_age_varying, compute_pcr_mfi_density

dir.create("paper", showWarnings = FALSE)

ANTIGENS <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")
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
R_FIXED  <- 8

COL_PCRN <- "#377eb8"
COL_PCRP <- "#e41a1c"
COL_NEG  <- "blue"
COL_POS  <- "red"

results_subdir_global <- sprintf("skew_normal_fixed_r%.0f_1mu_rho_age", R_FIXED)

LEG_BAR_NEG  <- "Observed PCR-"
LEG_BAR_POS  <- "Observed PCR+"
LEG_LINE_TOT <- "Fitted skew-normal total mixture"
LEG_LINE_POS <- "Fitted skew-normal mixture, PCR+ weighted"

make_row1_1mu_pcr_weighted <- function(ag, log_mfi, pcr_bin, post, age) {
  N <- length(log_mfi)
  
  mu_neg   <- median(post$mu_neg)
  mu_pos   <- median(post$mu_pos)
  sig_neg  <- median(post$sigma_neg)
  sig_pos  <- median(post$sigma_pos)
  skew_neg <- median(post$skew_neg)
  skew_pos <- median(post$skew_pos)
  
  pcr_prev <- mean(pcr_bin)
  
  S_med <- sapply(age, function(a)
    calculate_compartments_age_varying(a, median(post$lambda), median(post$alpha),
                                       apply(post$rho, 2, median), R_FIXED))
  
  bw     <- (max(log_mfi) - min(log_mfi)) / 50
  breaks <- seq(min(log_mfi) - bw / 2, max(log_mfi) + bw * 1.5, by = bw)
  
  h_neg <- hist(log_mfi[pcr_bin == 0], breaks = breaks, plot = FALSE)
  h_pos <- hist(log_mfi[pcr_bin == 1], breaks = breaks, plot = FALSE)
  
  df_hist <- rbind(
    data.frame(x = h_neg$mids, y = h_neg$counts / (N * bw), pcr_lab = LEG_BAR_NEG),
    data.frame(x = h_pos$mids, y = h_pos$counts / (N * bw), pcr_lab = LEG_BAR_POS)
  )
  df_hist$pcr_lab <- factor(df_hist$pcr_lab, levels = c(LEG_BAR_NEG, LEG_BAR_POS))
  
  xg   <- seq(min(log_mfi) - 0.5, max(log_mfi) + 0.5, length.out = 500)
  dens <- compute_pcr_mfi_density(S_med, pcr_bin, xg,
                                  mu_neg, sig_neg, skew_neg,
                                  mu_pos, sig_pos, skew_pos)
  
  d_pcrneg <- (1 - pcr_prev) * dens$d_pcrneg
  d_pcrpos <- pcr_prev       * dens$d_pcrpos
  d_total  <- d_pcrneg + d_pcrpos
  
  df_neg_line <- data.frame(x = xg, y = d_total, grp = LEG_LINE_TOT)  # PCR- + PCR+ (blue), from 0
  df_pos_line <- data.frame(x = xg, y = d_pcrpos, grp = LEG_LINE_POS) # PCR+ (red), stacked below of blue
  
  ggplot() +
    geom_col(data = df_hist,
             aes(x = x, y = y, fill = pcr_lab),
             width = bw * 0.85, position = "stack", alpha = 0.55) +
    scale_fill_manual(values = setNames(c(COL_PCRN, COL_PCRP), c(LEG_BAR_NEG, LEG_BAR_POS)),
                      name = NULL) +
    geom_line(data = df_neg_line, aes(x = x, y = y, colour = grp),
              linewidth = 1.1) +
    geom_line(data = df_pos_line, aes(x = x, y = y, colour = grp),
              linewidth = 1.2) +
    scale_colour_manual(name = NULL, values = setNames(c(COL_NEG, COL_POS), c(LEG_LINE_TOT, LEG_LINE_POS))) +
    labs(x = "log(MFI)", y = "Density", title = ag) +
    theme_classic(base_size = 26) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 30),
      axis.title = element_text(size = 26, face = "bold"),
      axis.text  = element_text(size = 26)
    )
}

panels <- vector("list", length(ANTIGENS))
for (i in seq_along(ANTIGENS)) {
  ag  <- ANTIGENS[i]
  mdr <- file.path("results", ag, results_subdir_global)
  ps <- readRDS(file.path(mdr, "posterior_samples.rds"))
  sd <- readRDS(file.path(mdr, "stan_data.rds"))
  panels[[i]] <- make_row1_1mu_pcr_weighted(
    ag      = ag,
    log_mfi = sd$log_mfi,
    pcr_bin = sd$pcr,
    post    = ps,
    age     = sd$age
  )
}

fig4 <- wrap_plots(panels, nrow = 2, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 24),
        legend.key.size = unit(1.4, "cm"),
        legend.box      = "vertical")

ggsave("paper/figure4_1mu_rho_age_hist_pcr_weighted.png", fig4,
       width = 30, height = 17, dpi = 300, bg = "white")
