data{
  int<lower=1> num_tracts;
  int<lower=1> num_vars;
  int<lower=1> num_years;
  matrix[num_years, num_vars] tract_m[num_tracts];
  vector[num_years] obs_resp[num_tracts];
  int resp_is_obs[num_tracts, num_years];
}
transformed data {
  vector[num_years] gp_mu;
  real years[num_years];
  gp_mu = rep_vector(0, num_years);
  for (i in 1:num_years)
    years[i] = i;
}
parameters{
  vector[num_vars] betas[num_years];
  real<lower = 0> rho[num_tracts];
  real<lower = 0> alpha[num_tracts];
  real<lower = 0> sigma[num_tracts];
}
model{
  vector[num_years] mu;
  matrix[num_years, num_years] K;
  matrix[num_years, num_years] L_K;
  real sq_sigma;
  
  for (i in 1:num_years) {
    betas[i] ~ normal(0, 10);
  }
  for (i in 1:num_tracts) {
    for (j in 1:num_years) {
      mu[j] = tract_m[i][j] * betas[j];
      if (resp_is_obs[i, j] == 1) {
        mu[j] ~ normal(obs_resp[i][j], 10);
      } 
      else if (resp_is_obs[i, j] == -1) {
        mu[j] = obs_resp[i][j];
      }
    }
    
    K = cov_exp_quad(years, alpha[i], rho[i]);
    sq_sigma = square(sigma[i]);
    for (k in 1:num_years) {
      K[k, k] = K[k, k] + sq_sigma;
    }
    L_K = cholesky_decompose(K);
    rho[i] ~ gamma(4, 4);
    alpha[i] ~ normal(0, 1);
    sigma[i] ~ normal(0, 1);
    
    mu ~ multi_normal_cholesky(gp_mu, L_K);
  }
}
generated quantities{
  vector[num_tracts] mu[num_years];
  for (i in 1:num_tracts) {
    for (j in 1:num_years) {
      mu[j][i] = tract_m[i][j] * betas[j];
    }
  }
}
