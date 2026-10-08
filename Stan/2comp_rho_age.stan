////////////////////////////////////////////////////////////////////////////////
// SEROCATALYTIC MODEL — 2 COMPARTMENTS
//
//   1. NO r (PCR clearance not modelled)
//   2. SKEW-NORMAL MFI distributions
//   3. AGE-VARYING rho — CHAINED PROPAGATORS (3 groups)
//        Group 1:  0–5  years  (rho[1])
//        Group 2:  5–18 years  (rho[2])
//        Group 3:  >18 years   (rho[3])
//   4. SINGLE mu_pos for seropositive class
//
//   ODE:
//     dS0/dt = -lambda*S0 + rho*S1
//     dS1/dt =  lambda*S0 - rho*S1
//
//   Explicit 2x2 propagator (exact closed-form solution):
//     Phi(tau, lambda, rho) = 1/(lambda+rho) *
//       [[rho + lambda*E,  rho*(1-E)  ],
//        [lambda*(1-E),    lambda + rho*E]]
//     where E = exp(-(lambda+rho)*tau)
//
//   All individuals contribute to the MFI likelihood (PCR not used).
//   Likelihood: log_sum_exp(log(S0) + f_neg(mfi), log(S1) + f_pos(mfi))
////////////////////////////////////////////////////////////////////////////////

functions {

  // Exact 2x2 propagator for the 2-compartment serocatalytic ODE
  vector propagate_2comp(real tau, real lambda, real rho, vector S_init) {
    real denom = lambda + rho;
    real E     = exp(-denom * tau);
    real phi00 = (rho + lambda * E) / denom;
    real phi01 = rho * (1.0 - E) / denom;
    real phi10 = lambda * (1.0 - E) / denom;
    real phi11 = (lambda + rho * E) / denom;
    vector[2] S;
    S[1] = phi00 * S_init[1] + phi01 * S_init[2];
    S[2] = phi10 * S_init[1] + phi11 * S_init[2];
    return S;
  }

  // Chained propagator — age boundaries mirror compartment_model.R: a1 and a2
  vector calculate_compartments_age_varying(
      real age,
      real lambda,
      vector rho
  ) {
    real a1 = 5.0;   // first age boundary (years)
    real a2 = 18.0;   // second age boundary (years)
    vector[2] S_birth = [1.0, 0.0]';

    if (age <= a1) {
      return propagate_2comp(age, lambda, rho[1], S_birth);

    } else if (age <= a2) {
      vector[2] S_at_a1 = propagate_2comp(a1, lambda, rho[1], S_birth);
      return propagate_2comp(age - a1, lambda, rho[2], S_at_a1);

    } else {
      vector[2] S_at_a1 = propagate_2comp(a1,      lambda, rho[1], S_birth);
      vector[2] S_at_a2 = propagate_2comp(a2 - a1, lambda, rho[2], S_at_a1);
      return propagate_2comp(age - a2, lambda, rho[3], S_at_a2);
    }
  }

}

// ----------------------------------------------------------------------------
data {
  int<lower=0> N;
  vector<lower=0>[N] age;
  vector[N] log_mfi;
  array[N] int<lower=0, upper=1> pcr;  // present for data consistency, not used in likelihood

  real<lower=0> lambda_max;
}

// ----------------------------------------------------------------------------
transformed data {
  real log_mfi_min   = min(log_mfi);
  real log_mfi_max   = max(log_mfi);
  real log_mfi_range = log_mfi_max - log_mfi_min;
}

// ----------------------------------------------------------------------------
parameters {
  real<lower=0, upper=2>      lambda;
  vector<lower=0, upper=6>[3] rho;

  real<lower=log_mfi_min, upper=log_mfi_max> mu_neg;
  real<lower=0>  sigma_neg;
  real           skew_neg;

  real<lower=4, upper=log_mfi_range> delta_sero;
  real<lower=0>  sigma_pos;
  real           skew_pos;
}

// ----------------------------------------------------------------------------
transformed parameters {
  real mu_pos = mu_neg + delta_sero;
}

// ----------------------------------------------------------------------------
model {
  lambda ~ exponential(10.0);
  rho    ~ exponential(10.0);

  mu_neg    ~ normal(4.0, 1.0);
  sigma_neg ~ exponential(1.0);
  skew_neg  ~ normal(0, 5);

  delta_sero ~ normal(5.0, 0.8);
  sigma_pos  ~ exponential(1.0);
  skew_pos   ~ normal(0, 5);

  for (i in 1:N) {
    vector[2] S = calculate_compartments_age_varying(age[i], lambda, rho);

    real prob_S0 = fmax(S[1], 1e-10);
    real prob_S1 = fmax(S[2], 1e-10);

    real ll_S0 = log(prob_S0) + skew_normal_lpdf(log_mfi[i] | mu_neg, sigma_neg, skew_neg);
    real ll_S1 = log(prob_S1) + skew_normal_lpdf(log_mfi[i] | mu_pos, sigma_pos, skew_pos);
    target += log_sum_exp(ll_S0, ll_S1);
  }
}

// ----------------------------------------------------------------------------
generated quantities {
  vector[N] prob_seropositive;
  vector[N] log_lik;

  for (i in 1:N) {
    vector[2] S = calculate_compartments_age_varying(age[i], lambda, rho);

    prob_seropositive[i] = S[2];

    real prob_S0 = fmax(S[1], 1e-10);
    real prob_S1 = fmax(S[2], 1e-10);

    real ll_S0 = log(prob_S0) + skew_normal_lpdf(log_mfi[i] | mu_neg, sigma_neg, skew_neg);
    real ll_S1 = log(prob_S1) + skew_normal_lpdf(log_mfi[i] | mu_pos, sigma_pos, skew_pos);
    log_lik[i] = log_sum_exp(ll_S0, ll_S1);
  }
}
