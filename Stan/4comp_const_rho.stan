////////////////////////////////////////////////////////////////////////////////
// SEROCATALYTIC-PCR MODEL WITH 4 COMPARTMENTS - FIXED R - SKEW-NORMAL MFI
////////////////////////////////////////////////////////////////////////////////
//
// ODE system:
//   dS0/dt = -lambda*S0 + r*S1 + rho*S3
//   dS1/dt = lambda*(1-alpha)*S0 - (alpha*lambda + r)*S1
//   dS2/dt = alpha*lambda*S0 + alpha*lambda*S1 - r*S2 + lambda*S3
//   dS3/dt = r*S2 - (rho + lambda)*S3
//
////////////////////////////////////////////////////////////////////////////////

functions {

  vector calculate_compartments(real t, real lambda, real alpha, real rho, real r) {
    vector[4] S;
    real mu2 = -(lambda + r);
    real b_coef = (1 + alpha) * lambda + r + rho;
    real c_coef = alpha * lambda * (lambda + r + rho) + r * rho;
    real Delta = b_coef * b_coef - 4 * c_coef;
    if (Delta < 0) Delta = 0;
    real sqrt_Delta = sqrt(Delta);
    real mu3 = (-b_coef + sqrt_Delta) / 2;
    real mu4 = (-b_coef - sqrt_Delta) / 2;
    real denom1 = alpha * lambda + r;
    real v21 = lambda * (1 - alpha) / denom1;
    real v41 = alpha * lambda * (lambda + r) / (rho * denom1);
    real v31 = (rho + lambda) * v41 / r;
    real v23 = lambda * (1 - alpha) / (mu3 + alpha * lambda + r);
    real numer3 = alpha * lambda * (1 + v23) * (mu3 + rho + lambda);
    real denom3 = (mu3 + r) * (mu3 + rho + lambda) - lambda * r;
    real v33 = numer3 / denom3;
    real v43 = r * v33 / (mu3 + rho + lambda);
    real v24 = lambda * (1 - alpha) / (mu4 + alpha * lambda + r);
    real numer4 = alpha * lambda * (1 + v24) * (mu4 + rho + lambda);
    real denom4 = (mu4 + r) * (mu4 + rho + lambda) - lambda * r;
    real v34 = numer4 / denom4;
    real v44 = r * v34 / (mu4 + rho + lambda);
    matrix[4, 4] V;
    V[1, 1] = 1;    V[1, 2] = 1;    V[1, 3] = 1;    V[1, 4] = 1;
    V[2, 1] = v21;  V[2, 2] = -1;   V[2, 3] = v23;  V[2, 4] = v24;
    V[3, 1] = v31;  V[3, 2] = 0;    V[3, 3] = v33;  V[3, 4] = v34;
    V[4, 1] = v41;  V[4, 2] = 0;    V[4, 3] = v43;  V[4, 4] = v44;
    vector[4] x0 = [1, 0, 0, 0]';
    vector[4] C = inverse(V) * x0;
    real exp_mu2_t = exp(mu2 * t);
    real exp_mu3_t = exp(mu3 * t);
    real exp_mu4_t = exp(mu4 * t);
    S[1] = C[1]*V[1,1] + C[2]*exp_mu2_t*V[1,2] + C[3]*exp_mu3_t*V[1,3] + C[4]*exp_mu4_t*V[1,4];
    S[2] = C[1]*V[2,1] + C[2]*exp_mu2_t*V[2,2] + C[3]*exp_mu3_t*V[2,3] + C[4]*exp_mu4_t*V[2,4];
    S[3] = C[1]*V[3,1] + C[2]*exp_mu2_t*V[3,2] + C[3]*exp_mu3_t*V[3,3] + C[4]*exp_mu4_t*V[3,4];
    S[4] = C[1]*V[4,1] + C[2]*exp_mu2_t*V[4,2] + C[3]*exp_mu3_t*V[4,3] + C[4]*exp_mu4_t*V[4,4];
    for (i in 1:4) if (S[i] < 0) S[i] = 0;
    real sum_S = S[1] + S[2] + S[3] + S[4];
    if (sum_S > 0) for (i in 1:4) S[i] = S[i] / sum_S;
    return S;
  }

  real calculate_seroprevalence(real t, real lambda, real alpha, real rho, real r) {
    vector[4] S = calculate_compartments(t, lambda, alpha, rho, r);
    return S[3] + S[4];
  }

  real calculate_pcr_prevalence(real t, real lambda, real alpha, real rho, real r) {
    vector[4] S = calculate_compartments(t, lambda, alpha, rho, r);
    return S[2] + S[3];
  }
}

data {
  int<lower=0> N;
  vector<lower=0>[N] age;
  vector[N] log_mfi;
  int<lower=0, upper=1> pcr[N];

  // r is FIXED as data (from literature)
  real<lower=0> r_fixed;

  // Bounds
  real<lower=0> lambda_max;
}

transformed data {
  real log_mfi_min = min(log_mfi);
  real log_mfi_max = max(log_mfi);
  real log_mfi_range = log_mfi_max - log_mfi_min;
  real r = r_fixed;
}

parameters {
  real<lower=0, upper=2> lambda;
  real<lower=0, upper=1> alpha;
  real<lower=0, upper=6> rho;

  // MFI distribution parameters - SKEW-NORMAL
  real<lower=log_mfi_min, upper=log_mfi_max> mu_neg;
  real<lower=4, upper=log_mfi_range> delta_mu;
  real<lower=0> sigma_neg;
  real<lower=0> sigma_pos;

  real skew_neg;
  real skew_pos;
}

transformed parameters {
  real mu_pos = mu_neg + delta_mu;

  vector[N] seroprevalence;
  vector[N] pcr_prevalence;

  for (i in 1:N) {
    seroprevalence[i] = calculate_seroprevalence(age[i], lambda, alpha, rho, r);
    pcr_prevalence[i] = calculate_pcr_prevalence(age[i], lambda, alpha, rho, r);

    if (seroprevalence[i] < 0) seroprevalence[i] = 0;
    if (seroprevalence[i] > 1) seroprevalence[i] = 1;
    if (pcr_prevalence[i] < 0) pcr_prevalence[i] = 0;
    if (pcr_prevalence[i] > 1) pcr_prevalence[i] = 1;
  }
}

model {
  // Priors
  lambda ~ exponential(10.0);
  alpha ~ beta(5, 2);
  rho ~ exponential(10.0);

  // MFI priors
  mu_neg ~ normal(4.0, 1);
  delta_mu ~ normal(5.0, 0.8);
  sigma_neg ~ exponential(1.0);
  sigma_pos ~ exponential(1.0);
  skew_neg ~ normal(0, 5);
  skew_pos ~ normal(0, 5);

  // Likelihood
  for (i in 1:N) {

    vector[4] S = calculate_compartments(age[i], lambda, alpha, rho, r);
    real prob_S0 = S[1];
    real prob_S1 = S[2];
    real prob_S2 = S[3];
    real prob_S3 = S[4];

    if (prob_S0 < 1e-10) prob_S0 = 1e-10;
    if (prob_S1 < 1e-10) prob_S1 = 1e-10;
    if (prob_S2 < 1e-10) prob_S2 = 1e-10;
    if (prob_S3 < 1e-10) prob_S3 = 1e-10;

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

generated quantities {
  vector[N] prob_seropositive;
  vector[N] prob_pcr_positive;
  vector[N] log_lik;

  real pcr_duration_days = 365.0 / r;

  for (i in 1:N) {
    vector[4] S = calculate_compartments(age[i], lambda, alpha, rho, r);

    prob_seropositive[i] = S[3] + S[4];
    prob_pcr_positive[i] = S[2] + S[3];

    real prob_S0 = S[1];
    real prob_S1 = S[2];
    real prob_S2 = S[3];
    real prob_S3 = S[4];

    if (prob_S0 < 1e-10) prob_S0 = 1e-10;
    if (prob_S1 < 1e-10) prob_S1 = 1e-10;
    if (prob_S2 < 1e-10) prob_S2 = 1e-10;
    if (prob_S3 < 1e-10) prob_S3 = 1e-10;

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
