################################################################################
# CANONICAL 2-COMPARTMENT MODEL 
#
# Stan source:
#   stan_models/serocatalytic_2compartments_skewnormal_rho_age.stan
#
# Compartment ordering:
#   S[1] = S0  : seronegative
#   S[2] = S1  : seropositive
#
# Seroprevalence = S1 = S["S1"]
# PCR not modelled in this version.
#
# Age-group boundaries are inherited from compartment_model.R (a1, a2).
#
# Source this file from the project root:
#   source("compartment_model_2comp.R")
################################################################################

source("compartment_model.R")   # provides a1, a2

# Exact 2x2 propagator: S(tau) = Phi(tau; lambda, rho) * S_init
propagate_2comp <- function(tau, lambda, rho, S_init) {
  S_init <- unname(S_init)   # prevent name-concatenation artefacts (e.g. "S0.S0")
  denom <- lambda + rho
  E     <- exp(-denom * tau)
  phi00 <- (rho + lambda * E) / denom
  phi01 <- rho * (1 - E) / denom
  phi10 <- lambda * (1 - E) / denom
  phi11 <- (lambda + rho * E) / denom
  c(
    S0 = phi00 * S_init[1] + phi01 * S_init[2],
    S1 = phi10 * S_init[1] + phi11 * S_init[2]
  )
}

# Chained propagator across age boundaries a1 and a2 (defined in compartment_model.R)
calculate_compartments_age_varying_2comp <- function(age, lambda, rho_vec) {
  S_birth <- c(S0 = 1, S1 = 0)

  if (age <= a1) {
    return(propagate_2comp(age,        lambda, rho_vec[1], S_birth))
  } else if (age <= a2) {
    S_at_a1 <- propagate_2comp(a1,     lambda, rho_vec[1], S_birth)
    return(propagate_2comp(age - a1,   lambda, rho_vec[2], S_at_a1))
  } else {
    S_at_a1 <- propagate_2comp(a1,      lambda, rho_vec[1], S_birth)
    S_at_a2 <- propagate_2comp(a2 - a1, lambda, rho_vec[2], S_at_a1)
    return(propagate_2comp(age - a2,    lambda, rho_vec[3], S_at_a2))
  }
}

# Convenience wrapper: from birth (S_init = (1,0)), constant rho
calculate_2comp <- function(t, lambda, rho) {
  propagate_2comp(t, lambda, rho, c(S0 = 1, S1 = 0))
}

################################################################################
# MFI UTILITY FUNCTIONS — require sn package (library(sn))
#
# S argument: 2 × N named matrix (rows S0, S1) produced by
#   sapply(ages, calculate_2comp(...)) or
#   sapply(ages, calculate_compartments_age_varying_2comp(...))
# pcr argument: integer vector length N  (0 = PCR−, 1 = PCR+)
#
# The 2-comp model has no PCR compartment.  PCR-stratified panels use the
# mean seroprevalence within each PCR subgroup (age-averaged marginal mix),
# not a PCR-conditioned probability as in the 4-comp model.
################################################################################

compute_2comp_mfi_density <- function(S, pcr, x_seq,
                                       mu_neg, sigma_neg, skew_neg,
                                       mu_pos, sigma_pos, skew_pos) {
  pi_vec   <- S["S1", ]
  w_pcrneg <- mean(pi_vec[pcr == 0])
  w_pcrpos <- mean(pi_vec[pcr == 1])
  pi_sero  <- mean(pi_vec)
  f_neg    <- dsn(x_seq, xi = mu_neg, omega = sigma_neg, alpha = skew_neg)
  f_pos    <- dsn(x_seq, xi = mu_pos, omega = sigma_pos, alpha = skew_pos)
  list(
    w_pcrneg = w_pcrneg,
    w_pcrpos = w_pcrpos,
    pi_sero  = pi_sero,
    d_pcrneg = (1 - w_pcrneg) * f_neg + w_pcrneg * f_pos,
    d_pcrpos = (1 - w_pcrpos) * f_neg + w_pcrpos * f_pos,
    d_all    = (1 - pi_sero)  * f_neg + pi_sero  * f_pos
  )
}

simulate_2comp_mfi <- function(S, mu_neg, sigma_neg, skew_neg,
                                mu_pos, sigma_pos, skew_pos) {
  pi_vec     <- S["S1", ]
  n          <- length(pi_vec)
  is_seropos <- rbinom(n, 1, pi_vec)
  ifelse(
    is_seropos == 1,
    rsn(n, xi = mu_pos, omega = sigma_pos, alpha = skew_pos),
    rsn(n, xi = mu_neg, omega = sigma_neg, alpha = skew_neg)
  )
}
