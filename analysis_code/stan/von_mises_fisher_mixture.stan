//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> D;
  int<lower=0> K;
  matrix[N, D] y;
  int<lower=0> N_pred;
  matrix[N_pred, D] y_pred;
  // matrix[K, D] prior_mean;
}

transformed data {
  int model_D = D;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  matrix[K, D] mu_tilde;
  vector<lower=0>[K] kappa;
  simplex[K] mixing_prob;
}

transformed parameters {
  matrix[K, D] mu;
  
  for (k in 1:K) {
    mu[k] = mu_tilde[k] / norm2(mu_tilde[k]);
  }
  
  // if (K == 1) {
  //   mixing_prob = 1;
  // }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  matrix[N, K] angle_prob = y * mu';
  
  if (K > 1) {
    mixing_prob ~ dirichlet(rep_vector(1, K));
  }
  
  for (k in 1:K) {
    mu_tilde[k] ~ normal(0, 1);
  }
  // for (k in 1:K) {
  //   real mu_tilde_norm = norm2(mu_tilde[k]);
  //   target += log_determinant(
  //     identity_matrix(D) / mu_tilde_norm -
  //     mu_tilde[k]' * mu_tilde[k] / pow(mu_tilde_norm, 3));
  // }
  kappa ~ gamma(1,1);
  for (n in 1:N) {
    vector[K] log_prob = rep_vector(0, K);
    if (K > 1) {
      log_prob = log(mixing_prob);
    }
    for (k in 1:K) {
      log_prob[k] += 
        angle_prob[n, k] * kappa[k] + 
        (model_D / 2.0 - 1) * log(kappa[k]) - 
        log_modified_bessel_first_kind(model_D / 2.0 - 1, kappa[k]) -
        model_D / 2.0 * log(2 * pi());
    }
    target += log_sum_exp(log_prob);
  }
}

generated quantities {
  matrix[N_pred, K] angle_pred_prob = y_pred * mu';
  real predict_lp = 0;
  
  for (n in 1:N_pred) {
    vector[K] log_prob = rep_vector(0, K);
    if (K > 1) {
      log_prob = log(mixing_prob);
    }
    for (k in 1:K) {
      log_prob[k] += 
        angle_pred_prob[n, k] * kappa[k] + 
        (model_D / 2.0 - 1) * log(kappa[k]) - 
        log_modified_bessel_first_kind(model_D / 2.0 - 1, kappa[k]) 
        - model_D / 2.0 * log(2 * pi());;
    }
    predict_lp += log_sum_exp(log_prob);
  }
}

