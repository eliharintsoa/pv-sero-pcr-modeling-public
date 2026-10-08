////////////////////////////////////////////////////////////////////////////////
// SEROCATALYTIC MODEL — 2 COMPARTMENTS — SKEW-NORMAL MFI — DISTRICT
//
// Classic two-state model: Seronegative <-> Seropositive
//   PM(t) = lambda/(lambda+rho) * (1 - exp(-(lambda+rho)*t))
//
// MFI mixture:
//   log(RAU) ~ PM * SN(mu_0+mu_1, epsilon_1, alpha_1)
//           + (1-PM) * SN(mu_0,   epsilon_0, alpha_0)
//
// 
////////////////////////////////////////////////////////////////////////////////

data {
  int<lower=0> N;
  array[N] real<lower=0> age;
  array[N] real          Y;       // log(RAU / MFI)
  real Ab_min;
  real Ab_max;
}

parameters {
  real<lower=0>             lambda;
  real<lower=0>             rho;

  real<lower=Ab_min, upper=Ab_max> mu_0;
  real<lower=4>               mu_1;      // delta: muPos = mu_0 + mu_1
  real<lower=0>                    epsilon_0;
  real<lower=0>                    epsilon_1;
  real                             alpha_0;
  real                             alpha_1;
}

transformed parameters {
  array[N] real PM;
  for (i in 1:N) {
    PM[i] = (lambda / (lambda + rho)) * (1 - exp(-(lambda + rho) * age[i]));
  }
}

model {
  real muPos = mu_0 + mu_1;
  array[2] real LPS;

  // Priors
  lambda    ~ exponential(10.0);
  rho       ~ exponential(10.0);
  mu_0      ~ normal(4, 1); #mean seroneg
  mu_1      ~ normal(5, 0.8);  #delta 
  epsilon_0 ~ exponential(1.0);
  epsilon_1 ~ exponential(1.0);
  alpha_0   ~ normal(0, 5);
  alpha_1   ~ normal(0, 5);

  // Likelihood
  for (j in 1:N) {
    real I1 = PM[j];
    real I0 = 1.0 - PM[j];

    LPS[1] = (I1 <= 0) ? -1e10 :
             log(I1) + skew_normal_lpdf(Y[j] | muPos, epsilon_1, alpha_1);
    LPS[2] = (I0 <= 0) ? -1e10 :
             log(I0) + skew_normal_lpdf(Y[j] | mu_0,  epsilon_0, alpha_0);

    target += log_sum_exp(LPS);
  }
}

generated quantities {
  array[N] real prob_seropositive;
  array[N] real log_lik;
  real muPos = mu_0 + mu_1;

  for (j in 1:N) {
    prob_seropositive[j] = PM[j];

    real I1 = PM[j];
    real I0 = 1.0 - PM[j];

    array[2] real LPS;
    LPS[1] = (I1 <= 0) ? -1e10 :
             log(I1) + skew_normal_lpdf(Y[j] | muPos, epsilon_1, alpha_1);
    LPS[2] = (I0 <= 0) ? -1e10 :
             log(I0) + skew_normal_lpdf(Y[j] | mu_0,  epsilon_0, alpha_0);

    log_lik[j] = log_sum_exp(LPS);
  }
}
