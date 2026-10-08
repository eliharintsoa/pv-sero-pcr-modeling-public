################################################################################
# CANONICAL 4-COMPARTMENT MODEL 
# Stan source:
#   stan_models/serocatalytic_pcr_4compartments_FIXED_R_skewnormal_commune.stan
#
# Compartment ordering (same as Stan S[1..4]):
#   S[1] = N0  : seroneg, PCR-
#   S[2] = N1  : seroneg, PCR+  (recent, pre-seroconversion)
#   S[3] = P1  : seropos, PCR+
#   S[4] = P0  : seropos, PCR-
#
# Seroprevalence = P1 + P0 = S["P1"] + S["P0"]
# PCR prevalence = N1 + P1 = S["N1"] + S["P1"]
#
# Age-group boundaries — change only here; all R scripts source this file.
#   Group 1: 0  – a1  years  →  rho[1]
#   Group 2: a1 – a2  years  →  rho[2]
#   Group 3: a2 – ∞   years  →  rho[3]
#
# Source this file from the project root:
#   source("compartment_model.R")
################################################################################

a1 <- 5    # first age boundary (years)
a2 <- 18    # second age boundary (years)

# Antigen display-name mapping — internal antigen code -> label shown in figures.
# Use ANTIGEN_LABELS[antigen_code] (or unname(...)) when labelling plots/facets,
# so every figure shows the same names without editing them one by one.
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

# General propagator: S(tau) = Phi(tau; lambda, alpha, rho, r) * S_init
# S_init: named or unnamed 4-vector (N0, N1, P1, P0)
calculate_compartments_from_init <- function(tau, lambda, alpha, rho, r, S_init) {

  mu2        <- -(lambda + r)
  b_coef     <- (1 + alpha) * lambda + r + rho
  c_coef     <- alpha * lambda * (lambda + r + rho) + r * rho
  Delta      <- b_coef * b_coef - 4 * c_coef
  if (Delta < 0) Delta <- 0
  sqrt_Delta <- sqrt(Delta)
  mu3 <- (-b_coef + sqrt_Delta) / 2
  mu4 <- (-b_coef - sqrt_Delta) / 2

  denom1 <- alpha * lambda + r
  v21 <- lambda * (1 - alpha) / denom1
  v41 <- alpha * lambda * (lambda + r) / (rho * denom1)
  v31 <- (rho + lambda) * v41 / r

  v23    <- lambda * (1 - alpha) / (mu3 + alpha * lambda + r)
  numer3 <- alpha * lambda * (1 + v23) * (mu3 + rho + lambda)
  denom3 <- (mu3 + r) * (mu3 + rho + lambda) - lambda * r
  v33    <- numer3 / denom3
  v43    <- r * v33 / (mu3 + rho + lambda)

  v24    <- lambda * (1 - alpha) / (mu4 + alpha * lambda + r)
  numer4 <- alpha * lambda * (1 + v24) * (mu4 + rho + lambda)
  denom4 <- (mu4 + r) * (mu4 + rho + lambda) - lambda * r
  v34    <- numer4 / denom4
  v44    <- r * v34 / (mu4 + rho + lambda)

  V <- matrix(c(
    1,   1,   1,   1,
    v21, -1,  v23, v24,
    v31,  0,  v33, v34,
    v41,  0,  v43, v44
  ), nrow = 4, byrow = TRUE)

  C <- solve(V, as.numeric(S_init))

  exp_mu2_t <- exp(mu2 * tau)
  exp_mu3_t <- exp(mu3 * tau)
  exp_mu4_t <- exp(mu4 * tau)

  S <- numeric(4)
  S[1] <- C[1]*V[1,1] + C[2]*exp_mu2_t*V[1,2] + C[3]*exp_mu3_t*V[1,3] + C[4]*exp_mu4_t*V[1,4]
  S[2] <- C[1]*V[2,1] + C[2]*exp_mu2_t*V[2,2] + C[3]*exp_mu3_t*V[2,3] + C[4]*exp_mu4_t*V[2,4]
  S[3] <- C[1]*V[3,1] + C[2]*exp_mu2_t*V[3,2] + C[3]*exp_mu3_t*V[3,3] + C[4]*exp_mu4_t*V[3,4]
  S[4] <- C[1]*V[4,1] + C[2]*exp_mu2_t*V[4,2] + C[3]*exp_mu3_t*V[4,3] + C[4]*exp_mu4_t*V[4,4]

  S[S < 0] <- 0
  sum_S <- sum(S)
  if (sum_S > 0) S <- S / sum_S

  c(N0 = S[1], N1 = S[2], P1 = S[3], P0 = S[4])
}

# Chained propagator across age boundaries a1 and a2 (defined at top of file)
calculate_compartments_age_varying <- function(age, lambda, alpha, rho_vec, r) {
  S_birth <- c(N0 = 1, N1 = 0, P1 = 0, P0 = 0)

  if (age <= a1) {
    return(calculate_compartments_from_init(age,        lambda, alpha, rho_vec[1], r, S_birth))
  } else if (age <= a2) {
    S_at_a1 <- calculate_compartments_from_init(a1,     lambda, alpha, rho_vec[1], r, S_birth)
    return(calculate_compartments_from_init(age - a1,   lambda, alpha, rho_vec[2], r, S_at_a1))
  } else {
    S_at_a1 <- calculate_compartments_from_init(a1,      lambda, alpha, rho_vec[1], r, S_birth)
    S_at_a2 <- calculate_compartments_from_init(a2 - a1, lambda, alpha, rho_vec[2], r, S_at_a1)
    return(calculate_compartments_from_init(age - a2,    lambda, alpha, rho_vec[3], r, S_at_a2))
  }
}

# Convenience wrapper: from birth (S_init = (1,0,0,0)), constant rho
calculate_compartments <- function(t, lambda, alpha, rho, r) {
  calculate_compartments_from_init(t, lambda, alpha, rho, r, c(1, 0, 0, 0))
}

################################################################################
# MFI UTILITY FUNCTIONS — require sn package (library(sn))
#
# S argument: 4 × N named matrix (rows N0, N1, P1, P0) produced by
#   sapply(ages, calculate_compartments(...)) or
#   sapply(ages, calculate_compartments_age_varying(...))
# pcr argument: integer vector length N  (0 = PCR−, 1 = PCR+)
################################################################################

# Analytical PCR-conditioned MFI density with individual-age-averaged weights.
# Weights are averaged over the actual age distribution of each PCR group 
#
# Returns a list:
#   w_pcrneg, w_pcrpos  — mean P(seropos | PCR group) over individuals
#   pi_sero             — overall mean seroprevalence
#   d_pcrneg, d_pcrpos, d_all — density vectors on x_seq
#   df                  — data.frame(x, density, pcr_status) ready for ggplot
compute_pcr_mfi_density <- function(S, pcr, x_seq,
                                    mu_neg, sigma_neg, skew_neg,
                                    mu_pos, sigma_pos, skew_pos) {
  pcrneg <- pcr == 0
  pcrpos <- pcr == 1

  w_pcrneg <- mean(S["P0", pcrneg] / pmax(S["N0", pcrneg] + S["P0", pcrneg], 1e-10))
  w_pcrpos <- mean(S["P1", pcrpos] / pmax(S["N1", pcrpos] + S["P1", pcrpos], 1e-10))
  pi_sero  <- mean(S["P1", ] + S["P0", ])

  dsn_neg  <- dsn(x_seq, xi = mu_neg, omega = sigma_neg, alpha = skew_neg)
  dsn_pos  <- dsn(x_seq, xi = mu_pos, omega = sigma_pos, alpha = skew_pos)

  d_pcrneg <- (1 - w_pcrneg) * dsn_neg + w_pcrneg * dsn_pos
  d_pcrpos <- (1 - w_pcrpos) * dsn_neg + w_pcrpos * dsn_pos
  d_all    <- (1 - pi_sero)  * dsn_neg + pi_sero  * dsn_pos

  list(
    w_pcrneg = w_pcrneg,
    w_pcrpos = w_pcrpos,
    pi_sero  = pi_sero,
    d_pcrneg = d_pcrneg,
    d_pcrpos = d_pcrpos,
    d_all    = d_all,
    df = data.frame(
      x          = rep(x_seq, 3),
      density    = c(d_pcrneg, d_pcrpos, d_all),
      pcr_status = factor(rep(c("PCR -", "PCR +", "All"), each = length(x_seq)),
                          levels = c("PCR -", "PCR +", "All"))
    )
  )
}

# Simulate MFI per individual: PCR-conditioned draw from seropos/seroneg mixture.
# For each person, p_sero = P(seropos | PCR status, age) from the compartment
# matrix, then rbinom → rsn from the appropriate component.
#
# Returns: numeric vector of length N (simulated log-MFI values)
simulate_mfi <- function(S, pcr, mu_neg, sigma_neg, skew_neg,
                                  mu_pos, sigma_pos, skew_pos) {
  N      <- ncol(S)
  p_sero <- numeric(N)
  pcrneg <- pcr == 0
  pcrpos <- pcr == 1

  p_sero[pcrneg] <- S["P0", pcrneg] / pmax(S["N0", pcrneg] + S["P0", pcrneg], 1e-10)
  p_sero[pcrpos] <- S["P1", pcrpos] / pmax(S["N1", pcrpos] + S["P1", pcrpos], 1e-10)

  is_seropos <- rbinom(N, 1, p_sero)
  ifelse(
    is_seropos == 1,
    rsn(N, xi = mu_pos, omega = sigma_pos, alpha = skew_pos),
    rsn(N, xi = mu_neg, omega = sigma_neg, alpha = skew_neg)
  )
}
