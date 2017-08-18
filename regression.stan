data{
  int<lower=1> num_obs;
  //int num_missing;
  int<lower=1> num_vars;
  matrix[num_obs, num_vars] obs_matrix;
  vector[num_obs] obs_resp;
  //matrix[num_missing, num_vars] missing_matrix;
}
parameters{
  vector[num_vars] betas;
  //vector[num_missing] missing_resp;
}
model{
  vector[num_obs] mu;
  //vector[num_missing] missing_mu;

  betas ~ normal(0, 10);
  mu = obs_matrix * betas;
  obs_resp ~ normal(mu, 10);
  //if (num_missing > 0) {
    //missing_mu = missing_matrix * betas;
    //missing_resp ~ normal(missing_mu, 10);
  //}
}
generated quantities{
  vector[num_obs] mu;
  mu = obs_matrix * betas;
}
