################################################################################
# MFI BOXPLOT BY PCR SPECIES STATUS — observed log(MFI), faceted by antigen
#
# x-axis: PCR status (Negative, Pf, Pv, Mixed)
# facets: one per antigen (EBP, LF005, LF010, LF016, MSP8, PTEX150, PvCSS, RBP2b.P87)
#
# Output: paper/mfi_boxplot_by_pcr_species.png / .pdf
################################################################################



suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

dir.create("paper", showWarnings = FALSE)

ANTIGENS <- c("EBP", "LF005", "LF010", "LF016", "MSP8", "PTEX150", "PvCSS", "RBP2b.P87")

PCR_LEVELS <- c("Negative", "Pf", "Pv", "Mixed")
PCR_COLOURS <- c(Negative = "grey60", Pf = "#FFDFBF", Pv = "#e41a1c", Mixed = "#377eb8")

# 1. Load and reshape data
data_clean <- read.csv("Data/data_cleaned.csv") |>
  filter(!is.na(pcr)) |>
  mutate(
    pcr_species = recode(pcr,
                         "Negative"             = "Negative",
                         "Pf infection"          = "Pf",
                         "Pv infection"          = "Pv",
                         "Pf & Pv co-infection"  = "Mixed"
    ),
    pcr_species = factor(pcr_species, levels = PCR_LEVELS)
  )

antigen_cols <- paste0("log_mfi_", ANTIGENS)
antigen_cols <- antigen_cols[antigen_cols %in% colnames(data_clean)]

df_long <- data_clean |>
  select(sample_id, pcr_species, all_of(antigen_cols)) |>
  pivot_longer(cols = all_of(antigen_cols),
               names_to = "Antigen", values_to = "log_MFI") |>
  mutate(
    Antigen = gsub("^log_mfi_", "", Antigen),
    Antigen = factor(Antigen, levels = ANTIGENS)
  ) |>
  filter(!is.na(log_MFI))

# 2. Plot
p <- ggplot(df_long, aes(x = pcr_species, y = log_MFI, fill = pcr_species)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4, linewidth = 0.3) +
  facet_wrap(~ Antigen, nrow = 2) +
  scale_fill_manual(values = PCR_COLOURS, name = "PCR status") +
  labs(
    title = "Observed log(MFI) distribution by PCR species status and antigen",
    x     = "PCR status",
    y     = "log(MFI)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 12),
    strip.text        = element_text(face = "bold", size = 9),
    strip.background  = element_rect(fill = "grey95", colour = NA),
    axis.text.x       = element_text(angle = 30, hjust = 1, size = 7),
    axis.text.y       = element_text(size = 7),
    axis.title        = element_text(size = 9),
    legend.position   = "none",
    panel.spacing.x   = unit(0.6, "lines"),
    panel.spacing.y   = unit(0.8, "lines")
  )

ggsave("paper/mfi_boxplot_by_pcr_species.png", p,
       width = 14, height = 7, dpi = 300, bg = "white")
ggsave("paper/mfi_boxplot_by_pcr_species.pdf", p,
       width = 14, height = 7)
cat("Saved: paper/mfi_boxplot_by_pcr_species.png\n")
cat("Saved: paper/mfi_boxplot_by_pcr_species.pdf\n")
