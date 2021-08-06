data {
  int<lower=0> N;
  vector[N] observed_dir;
  vector[N] trans_m_list;
}

parameters {
  real pi_param_mu;
  vector[N] pi_param_tilde_list;
}

model {
  vector[N] pi_param_list = exp(0.05 * pi_param_tilde_list + pi_param_mu);
  
  pi_param_mu ~ normal(0, 5);
  pi_param_tilde_list ~ normal(0, 1);
  for (n in 1:N) {
    if (pi_param_list[n] < 100) {
      observed_dir[n] ~ von_mises(trans_m_list[n], pi_param_list[n]);
    } else {
      observed_dir[n] ~ normal(trans_m_list[n], sqrt(1 / pi_param_list[n]));
    }
  }
}

