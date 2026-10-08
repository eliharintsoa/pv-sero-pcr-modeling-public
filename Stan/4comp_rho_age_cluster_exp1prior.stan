////////////////////////////////////////////////////////////////////////////////
// SEROCATALYTIC-PCR MODEL — 4 COMPARTMENTS
//
//   1. FIXED R
//   2. SKEW-NORMAL MFI distributions
//   3. AGE-VARYING rho — 3 groups
//        Group 1:  0–5  years  (rho[1])
//        Group 2:  5–18 years  (rho[2])
//        Group 3:  >18 years   (rho[3])
//   4. SINGLE mu_pos — shared MFI mean for S2 and S3
//        mu_pos = mu_neg + delta_sero
//   5. CLUSTER STRATIFICATION — one lambda per cluster
//        All other parameters shared across clusters
//
//   PRIOR VARIANT: lambda ~ exponential(1.0)
////////////////////////////////////////////////////////////////////////////////

functions {

  vector calculate_compartments_from_init(
      real tau,
      real lambda,
      real alpha,
      real rho,
      real r,
      vector S_init
  ) {
    vector[4] S;

    real mu2       = -(lambda + r);
    real b_coef    = (1 + alpha) * lambda + r + rho;
    real c_coef    = alpha * lambda * (lambda + r + rho) + r * rho;
    real Delta     = b_coef * b_coef - 4 * c_coef;
    if (Delta < 0) Delta = 0;
    real sqrt_Delta = sqrt(Delta);
    real mu3 = (-b_coef + sqrt_Delta) / 2;
    real mu4 = (-b_coef - sqrt_Delta) / 2;

    real denom1 = alpha * lambda + r;
    real v21 = lambda * (1 - alpha) / denom1;
    real v41 = alpha * lambda * (lambda + r) / (rho * denom1);
    real v31 = (rho + lambda) * v41 / r;

    real v23    = lambda * (1 - alpha) / (mu3 + alpha * lambda + r);
    real numer3 = alpha * lambda * (1 + v23) * (mu3 + rho + lambda);
    real denom3 = (mu3 + r) * (mu3 + rho + lambda) - lambda * r;
    real v33    = numer3 / denom3;
    real v43    = r * v33 / (mu3 + rho + lambda);

    real v24    = lambda * (1 - alpha) / (mu4 + alpha * lambda + r);
    real numer4 = alpha * lambda * (1 + v24) * (mu4 + rho + lambda);
    real denom4 = (mu4 + r) * (mu4 + rho + lambda) - lambda * r;
    real v34    = numer4 / denom4;
    real v44    = r * v34 / (mu4 + rho + lambda);

    matrix[4, 4] V;
    V[1, 1] = 1;    V[1, 2] = 1;    V[1, 3] = 1;    V[1, 4] = 1;
    V[2, 1] = v21;  V[2, 2] = -1;   V[2, 3] = v23;  V[2, 4] = v24;
    V[3, 1] = v31;  V[3, 2] = 0;    V[3, 3] = v33;  V[3, 4] = v34;
    V[4, 1] = v41;  V[4, 2] = 0;    V[4, 3] = v43;  V[4, 4] = v44;

    vector[4] C = inverse(V) * S_init;

    real exp_mu2_t = exp(mu2 * tau);
    real exp_mu3_t = exp(mu3 * tau);
    real exp_mu4_t = exp(mu4 * tau);

    S[1] = C[1]*V[1,1] + C[2]*exp_mu2_t*V[1,2] + C[3]*exp_mu3_t*V[1,3] + C[4]*exp_mu4_t*V[1,4];
    S[2] = C[1]*V[2,1] + C[2]*exp_mu2_t*V[2,2] + C[3]*exp_mu3_t*V[2,3] + C[4]*exp_mu4_t*V[2,4];
    S[3] = C[1]*V[3,1] + C[2]*exp_mu2_t*V[3,2] + C[3]*exp_mu3_t*V[3,3] + C[4]*exp_mu4_t*V[3,4];
    S[4] = C[1]*V[4,1] + C[2]*exp_mu2_t*V[4,2] + C[3]*exp_mu3_t*V[4,3] + C[4]*exp_mu4_t*V[4,4];

    for (i in 1:4) if (S[i] < 0) S[i] = 0;
    real sum_S = S[1] + S[2] + S[3] + S[4];
    if (sum_S > 0) for (i in 1:4) S[i] = S[i] / sum_S;

    return S;
  }

  // Age boundaries: a1 = 5, a2 = 18
  vector calculate_compartments_age_varying(
      real age,
      real lambda,
      real alpha,
      vector rho,
      real r
  ) {
    real a1 = 5.0;
    real a2 = 18.0;
    vector[4] S_birth = [1.0, 0.0, 0.0, 0.0]';

    if (age <= a1) {
      return calculate_compartments_from_init(
          age, lambda, alpha, rho[1], r, S_birth);

    } else if (age <= a2) {
      vector[4] S_at_a1 = calculate_compartments_from_init(
          a1, lambda, alpha, rho[1], r, S_birth);
      return calculate_compartments_from_init(
          age - a1, lambda, alpha, rho[2], r, S_at_a1);

    } else {
      vector[4] S_at_a1 = calculate_compartments_from_init(
          a1, lambda, alpha, rho[1], r, S_birth);
      vector[4] S_at_a2 = calculate_compartments_from_init(
          a2 - a1, lambda, alpha, rho[2], r, S_at_a1);
      return calculate_compartments_from_init(
          age - a2, lambda, alpha, rho[3], r, S_at_a2);
    }
  }

}

// ----------------------------------------------------------------------------
data {
  int<lower=0> N;
  vector<lower=0>[N] age;
  vector[N] log_mfi;
  array[N] int<lower=0, upper=1> pcr;

  // Cluster stratification
  int<lower=1> N_cluster;
  array[N] int<lower=1, upper=N_cluster> cluster;

  real<lower=0> r_fixed;
  real<lower=0> lambda_max;
}

// ----------------------------------------------------------------------------
transformed data {
  real log_mfi_min   = min(log_mfi);
  real log_mfi_max   = max(log_mfi);
  real log_mfi_range = log_mfi_max - log_mfi_min;
  real r             = r_fixed;
}

// ----------------------------------------------------------------------------
parameters {
  // One lambda per cluster; all other parameters shared
  vector<lower=0, upper=2>[N_cluster] lambda;
  real<lower=0, upper=1>      alpha;
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
  lambda     ~ exponential(1.0);
  alpha      ~ beta(5, 2);
  rho        ~ exponential(10.0);

  mu_neg     ~ normal(4.0, 1.0);
  sigma_neg  ~ exponential(1.0);
  skew_neg   ~ normal(0, 5);

  delta_sero ~ normal(5.0, 2.0);
  sigma_pos  ~ exponential(1.0);
  skew_pos   ~ normal(0, 5);

  for (i in 1:N) {
    vector[4] S = calculate_compartments_age_varying(
        age[i], lambda[cluster[i]], alpha, rho, r);

    real prob_S0 = fmax(S[1], 1e-10);
    real prob_S1 = fmax(S[2], 1e-10);
    real prob_S2 = fmax(S[3], 1e-10);
    real prob_S3 = fmax(S[4], 1e-10);

    if (pcr[i] == 0) {
      real ll_S0 = log(prob_S0) + skew_normal_lpdf(log_mfi[i] | mu_neg, sigma_neg, skew_neg);
      real ll_S3 = log(prob_S3) + skew_normal_lpdf(log_mfi[i] | mu_pos, sigma_pos, skew_pos);
      target += log_sum_exp(ll_S0, ll_S3);
    } else {
      real ll_S1 = log(prob_S1) + skew_normal_lpdf(log_mfi[i] | mu_neg, sigma_neg, skew_neg);
      real ll_S2 = log(prob_S2) + skew_normal_lpdf(log_mfi[i] | mu_pos, sigma_pos, skew_pos);
      target += log_sum_exp(ll_S1, ll_S2);
    }
  }
}

// ----------------------------------------------------------------------------
generated quantities {
  vector[N] prob_seropositive;
  vector[N] prob_pcr_positive;
  vector[N] log_lik;

  real pcr_duration_days = 365.0 / r_fixed;

  for (i in 1:N) {
    vector[4] S = calculate_compartments_age_varying(
        age[i], lambda[cluster[i]], alpha, rho, r);

    prob_seropositive[i] = S[3] + S[4];
    prob_pcr_positive[i] = S[2] + S[3];

    real prob_S0 = fmax(S[1], 1e-10);
    real prob_S1 = fmax(S[2], 1e-10);
    real prob_S2 = fmax(S[3], 1e-10);
    real prob_S3 = fmax(S[4], 1e-10);

    if (pcr[i] == 0) {
      real ll_S0 = log(prob_S0) + skew_normal_lpdf(log_mfi[i] | mu_neg, sigma_neg, skew_neg);
      real ll_S3 = log(prob_S3) + skew_normal_lpdf(log_mfi[i] | mu_pos, sigma_pos, skew_pos);
      log_lik[i] = log_sum_exp(ll_S0, ll_S3);
    } else {
      real ll_S1 = log(prob_S1) + skew_normal_lpdf(log_mfi[i] | mu_neg, sigma_neg, skew_neg);
      real ll_S2 = log(prob_S2) + skew_normal_lpdf(log_mfi[i] | mu_pos, sigma_pos, skew_pos);
      log_lik[i] = log_sum_exp(ll_S1, ll_S2);
    }
  }
}
