data {
  int<lower=1> N;                  // number of observations
  vector<lower=0,upper=1>[N] D;    // Dissimilarity (as response)
  vector[N] delta_date;            // date_index
  vector[N] time_cat;              // time_index (day | night)
  vector[N] tre_cat;               // treatment_index (light | no light)
  vector[N] rep_cat;               // biological_replicate_index
  int<lower=1> J;
  array[N] int beta_idx;               // biological_replicate_index
}

parameters {
  real gamma;                      // coefficient for delta_date
  real alpha;                      // coefficient for time_cat
  // real beta;                       // coefficient for tre_cat
  vector[J] beta;                       // coefficient for tre_cat
  real theta;                      // coefficient for rep_cat
  real intercept;                  // intercept
  real<lower=0> phi;               // precision (dispersion) parameter
}

transformed parameters {
  vector<lower=0,upper=1>[N] mu;
  
  mu = inv_logit(intercept
                 + gamma * delta_date
                 + alpha * time_cat
                 + beta[beta_idx] .* tre_cat
                 + theta * rep_cat);
}

model {
  // priors (can be adjusted)
  gamma ~ normal(0, 1);
  alpha ~ normal(0, 1);
  beta ~  normal(0, 1);
  theta ~ normal(0, 1);
  intercept ~ normal(0, 1);
  phi ~ exponential(1);

  D ~ beta(mu * phi, (1 - mu) * phi);
}

// generated quantities {
//   // vector[N] log_lik;             // pointwise log-likelihood
//   vector<lower=0,upper=1>[N] D_est;  // posterior predictive simulations
// 
//   for (n in 1:N) {
//     // log-likelihood of each observation
//     log_lik[n] = beta_lpdf(D[n] | mu[n] * phi, (1 - mu[n]) * phi);
// 
//     // posterior predictive sample
//     D_est[n] = beta_rng(mu[n] * phi, (1 - mu[n]) * phi);
//   }
// }