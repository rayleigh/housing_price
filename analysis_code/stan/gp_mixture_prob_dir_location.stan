functions {
  matrix kronecker_prod(matrix A, matrix B) {
    matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
    int m = rows(A);
    int n = cols(A);
    int p = rows(B);
    int q = cols(B);
    for (i in 1:m) {
      for (j in 1:n) {
        int row_start = (i - 1) * p + 1;
        int row_end = (i - 1) * p + p;
        int col_start = (j - 1) * q + 1;
        int col_end = (j - 1) * q + q;
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
    return C;
  }
  
  vector angle_mod(vector A, int v_size) {
    vector[v_size] mod_A = A;
    for (i in 1:v_size) {
      if (A[i] < 0) {
        mod_A[i] = A[i] + 2 * pi();
      } else if (A[i] > 2 * pi()) {
        mod_A[i] = A[i] - 2 * pi(); 
      }
    }
    return(mod_A);
  }
}

data {
  int<lower=0> N;
  int<lower=0> N_pred;
  int<lower=0> d;
  vector[N] observed_dir;
  row_vector[d] simplex_location_m[N];
  real<lower = 0> m_l;
  real<lower = 0> m_s;
  row_vector[d] simplex_prediction_m[N_pred];
  vector[N_pred] pred_dir;
  // real<lower = 0> pi_l;
}


transformed data {
  matrix[N, N] m_L_K;
  matrix[N_pred, N_pred] pred_cov_L_K;
  matrix[N, N_pred] k_m1_m2 = cov_exp_quad(simplex_location_m, simplex_prediction_m, m_s, m_l);
  real delta = 1e-9;


  {  
    matrix[N, N_pred] m_pre_cov_m;
    matrix[N_pred, N_pred] pred_cov;
    matrix[N, N] m_K = cov_exp_quad(simplex_location_m, m_s, m_l);
    matrix[N_pred, N_pred] pred_m_K = cov_exp_quad(simplex_prediction_m, m_s, m_l);
  
    for (n in 1:N) {
      m_K[n, n] += delta;
    }
    m_L_K = cholesky_decompose(m_K);
    
    m_pre_cov_m = mdivide_left_tri_low(m_L_K, k_m1_m2);
    pred_cov = pred_m_K - 
      m_pre_cov_m' * m_pre_cov_m + diag_matrix(rep_vector(delta, N_pred));
    pred_cov_L_K = cholesky_decompose(pred_cov);
  }
  
}

parameters {
  vector[N] p_tilde_list;
  ordered[2] ordered_logit_m_param_list;
  vector<lower = 0>[2] pi_param_list;
  //vector<lower = -2 * pi(), upper = 2 * pi()>[2] m_tilde_param_list;
}

transformed parameters {
  vector<lower = 0, upper = 1>[N] trans_p_list;
  vector<lower = -pi(), upper = 3 * pi()>[2] m_tilde_param_list = 4 * pi() * inv_logit(ordered_logit_m_param_list) - pi();
  //vector<lower = 0, upper = 2 * pi()>[2] m_param_list = 2 * pi() * inv_logit(ordered_logit_m_param_list);
  vector<lower = 0, upper = 2 * pi()>[2] m_param_list;
  
  {
    vector[N] p_list = m_L_K * p_tilde_list;
    trans_p_list = inv_logit(p_list);
  }
  m_param_list = angle_mod(m_tilde_param_list, 2);
}

model {
  // m_param_list ~ von_mises(pi(), 0.01);
  // m_param_list ~ uniform(0, 2 * pi());
  pi_param_list ~ gamma(1, 1);
  p_tilde_list ~ normal(0, 1);
  target += sum(ordered_logit_m_param_list - 2 * log1p_exp(ordered_logit_m_param_list));
  // target += 2 * log(2 * pi());

  for (n in 1:N) {
    if (observed_dir[n] == 0) {
      continue;
    }
    if (pi_param_list[1] < 100 && pi_param_list[2] < 100) {
      target += log_mix(trans_p_list[n],
                  von_mises_lpdf(observed_dir[n] | m_param_list[1], pi_param_list[1]),
                  von_mises_lpdf(observed_dir[n] | m_param_list[2], pi_param_list[2]));
    } else if (pi_param_list[1] >= 100 && pi_param_list[2] >= 100) {
      target += log_mix(trans_p_list[n],
                  normal_lpdf(observed_dir[n] | m_param_list[1], sqrt(1 / pi_param_list[1])),
                  normal_lpdf(observed_dir[n] | m_param_list[2], sqrt(1 / pi_param_list[2])));
    } else if (pi_param_list[1] >= 100) {
      target += log_mix(trans_p_list[n],
                  normal_lpdf(observed_dir[n] | m_param_list[1], sqrt(1 / pi_param_list[1])),
                  von_mises_lpdf(observed_dir[n] | m_param_list[2], pi_param_list[2]));
    } else {
      target += log_mix(trans_p_list[n],
                  von_mises_lpdf(observed_dir[n] | m_param_list[1], pi_param_list[1]),
                  normal_lpdf(observed_dir[n] | m_param_list[2], sqrt(1 / pi_param_list[2])));
    }
  }
}

generated quantities {
  real predict_lp[100];

  {
    real predict_lp_tmp;
    vector[N_pred] p_predict_tmp;
    vector[N] L_K_div_m_tilde = mdivide_right_tri_low(p_tilde_list', m_L_K)';
    vector[N_pred] pred_mu = k_m1_m2' * L_K_div_m_tilde;

    for (i in 1:size(predict_lp)) {
      predict_lp_tmp = 0;
      p_predict_tmp = inv_logit(multi_normal_cholesky_rng(pred_mu, pred_cov_L_K));
      for (n in 1:N_pred) {
        if (pi_param_list[1] < 100 && pi_param_list[2] < 100) {
          predict_lp_tmp += log_mix(p_predict_tmp[n],
                      von_mises_lpdf(pred_dir[n] | m_param_list[1], pi_param_list[1]),
                      von_mises_lpdf(pred_dir[n] | m_param_list[2], pi_param_list[2]));
        } else if (pi_param_list[1] >= 100 && pi_param_list[2] >= 100) {
          predict_lp_tmp += log_mix(p_predict_tmp[n],
                      normal_lpdf(pred_dir[n] | m_param_list[1], sqrt(1 / pi_param_list[1])),
                      normal_lpdf(pred_dir[n] | m_param_list[2], sqrt(1 / pi_param_list[2])));
        } else if (pi_param_list[1] >= 100) {
          predict_lp_tmp += log_mix(p_predict_tmp[n],
                      normal_lpdf(pred_dir[n] | m_param_list[1], sqrt(1 / pi_param_list[1])),
                      von_mises_lpdf(pred_dir[n] | m_param_list[2], pi_param_list[2]));
        } else {
          predict_lp_tmp += log_mix(p_predict_tmp[n],
                      von_mises_lpdf(pred_dir[n] | m_param_list[1], pi_param_list[1]),
                      normal_lpdf(pred_dir[n] | m_param_list[2], sqrt(1 / pi_param_list[2])));
        }
      }
      predict_lp[i] = predict_lp_tmp;
    }
  }
}
