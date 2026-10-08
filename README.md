# Pv_sero_pcr_modeling

Analysis code for a joint serocatalytic + PCR-prevalence Bayesian model of
*Plasmodium vivax* transmission, fitted to age-stratified antibody (MFI) and
PCR data from a community survey (Madagascar). Models are implemented in
Stan and fitted via `rstan`.

> **Status:** this repository is being populated incrementally. The model
> definitions, all 8 model-fitting scripts, all figure/table generation
> scripts (Figures 4–6, Table 2, r-sensitivity, plus 7 supplementary
> scripts), and `Data/` are complete.
> Table 1 (participant counts), Table 3, and the Ethiopia validation
> analyses are not included.

## Requirements

R (≥ 4.2) with:

```r
install.packages(c(
  "rstan", "bayesplot", "loo", "posterior",
  "ggplot2", "dplyr", "tidyr", "purrr", "tibble",
  "sn", "patchwork"
))
```

Scripts `01`–`11` assume the working directory is the repository root
(`Pv_sero_pcr_modeling/`) when run, e.g. `Rscript Script/02_run_stan_2comp_const_rho.R`
from that directory, or opening the folder as the RStudio project root.
Scripts `S1`–`S10` locate the repository root automatically from their own
file path, so they can be run (via `Rscript` or RStudio's Source button)
from any working directory, with no manual `setwd()` needed.

## Repository structure

```
Pv_sero_pcr_modeling/
├── compartment_model.R       # Shared ODE / density helper functions (4-comp + 2-comp)
├── compartment_model_2comp.R # 2-comp-specific helpers (sources compartment_model.R)
├── Stan/                     # Stan model definitions (.stan)
├── Data/                     # data_cleaned.csv (individual-level survey data)
└── Script/
    ├── 01_prepareStanData.R   # Builds Data/stan_data_*.rds from data_cleaned.csv
    ├── 02–06_runStan_*.R       # Model fitting (one script per model variant)
    ...
```

Run `Script/01_prepareStanData.R` once, first — it reads `Data/data_cleaned.csv`
and writes the `Data/stan_data_<antigen>.rds` and `Data/stan_data_<antigen>_stratified.rds`
files that every fitting script (02–06) and the r-sensitivity script (11) load via
`readRDS()`. Only `data_cleaned.csv` is checked into `Data/`; the `.rds` files are
generated, not stored, to avoid shipping 16 files that are fully derivable from
the one CSV.

## Model-fitting scripts

Each script fits **one antigen** per run. Set `antigen_name` before sourcing,
or edit the default at the top of the script, then run once per antigen
(`EBP`, `LF005`, `LF010`, `LF016`, `MSP8`, `PTEX150`, `PvCSS`, `RBP2b.P87`).
To fit all 8 antigens, wrap the script in a loop over `antigen_name` from the
outside, e.g.:

```r
for (ag in c("EBP","LF005","LF010","LF016","MSP8","PTEX150","PvCSS","RBP2b.P87")) {
  antigen_name <- ag
  source("Script/02_run_stan_2comp_const_rho.R")
}
```


