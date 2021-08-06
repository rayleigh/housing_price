library(LaplacesDemon)

angle <- function(x,y) {
  rotation = atan(y / x)
  if (x < 0) {
    return(rotation + pi)
  }
  if (y < 0) {
    return(rotation + 2 * pi)
  }
  return(rotation)
}

calc_von_mises_prob <- function(y, m_k, kappa_k, log = F) {
  vm_log_prob = kappa_k * cos(y - m_k) - log(2 * pi * besselI(kappa_k, 0))
  if (log) {
    return(vm_log_prob)
  }
  return(exp(vm_log_prob))
}

sample_m_i <- function(m_i, y_i, rho_i, r_i, mu_i1, mu_i2, sigma_i) {
  next_m_i = m_i + rnorm(1, 0, 0.01)
  if (next_m_i < 0) {
    next_m_i = next_m_i + 2 * pi
  } else if (next_m_i > 2 * pi) {
    next_m_i = next_m_i - 2 * pi
  }
  log_accept_prob = 
    calc_von_mises_prob(y_i, next_m_i, rho_i, log = T) +
    dnorm(r_i * cos(next_m_i), mu_i1, sigma_i, log = T) +
    dnorm(r_i * sin(next_m_i), mu_i2, sigma_i, log = T) - 
    (calc_von_mises_prob(y_i, m_i, rho_i, log = T) +
      dnorm(r_i * cos(m_i), mu_i1, sigma_i, log = T) +
      dnorm(r_i * sin(m_i), mu_i2, sigma_i, log = T)) 
  if (log(runif(1)) < log_accept_prob) {
    return(next_m_i)
  }
  return(m_i)
}

sample_r_i <- function(r_i, m_i, mu_i1, mu_i2, sigma_i) {
  next_r_i = r_i + rnorm(1, 0, 0.1)
  alpha_i = angle(mu_i1, mu_i2)
  mu_i0 = sqrt(mu_i1^2 + mu_i2^2)
  if (next_r_i < 0) {
    return(r_i)
  }
  log_accept_prob = 
    log(next_r_i) - 1 / (2 * sigma_i^2) * (next_r_i - mu_i0 * cos(m_i - alpha_i))^2 -
    (log(r_i) - 1 / (2 * sigma_i^2) * (r_i - mu_i0 * cos(m_i - alpha_i))^2)
  if (log(runif(1)) < log_accept_prob) {
    return(next_r_i)
  }
  return(r_i)
}

sample_z_tilde_gp_dir_slice_sampler <- function(prev_z_tilde_list, y, rho_list, gp_chol, mu_list) {
  nu <- matrix(rnorm(length(y) * ncol(prev_z_tilde_list)), ncol = 2)
  nu <- gp_chol %*% nu
  prev_m_list <- apply(prev_z_tilde_list, 1, function(row) {
    angle(row[1] + mu_list[1], row[2] + mu_list[2])
  })
  log_y <- sum(calc_von_mises_prob(y, prev_m_list, rho_list, log = T)) + log(runif(1))
  theta <- runif(1, 0, 2 * pi)
  theta_b <- c(theta - 2 * pi, theta)
  
  repeat {
    next_z_tilde_list <- prev_z_tilde_list * cos(theta) + nu * sin(theta)
    next_m_list <- apply(next_z_tilde_list, 1, function(row) {
      angle(row[1] + mu_list[1], row[2] + mu_list[2])
    })
    next_ll <- sum(calc_von_mises_prob(y, next_m_list, rho_list, log = T))
    if (next_ll > log_y) {
      break
    }
    if (theta < 0) {
      theta_b[1] = theta
    } else {
      theta_b[2] = theta
    }
    theta <- runif(1, theta_b[1], theta_b[2])
  }
  return(next_z_tilde_list)
}

calc_svm_c_ll <- function(y, p_list, z_tilde_list, rho_list, mu_list, K, comp_m_list = NULL) {
  y_tmp = rep(0, length(y))
  for (k in 1:K) {
    if (is.null(comp_m_list)) {
      m_list = mapply(function(tilde_z1, tilde_z2, mu_1, mu_2) {
        angle(tilde_z1 + mu_1, tilde_z2 + mu_2)
      }, 
      z_tilde_list[, 2 * (k - 1) + 1], z_tilde_list[, 2 * (k - 1) + 2],
      mu_list[2 * (k - 1) + 1], mu_list[2 * (k - 1) + 1])
    } else {
      m_list <- comp_m_list[, k]
    }
    y_tmp = log_y_tmp + p_list[k] * calc_von_mises_prob(y, m_list, rho_list[, k])
  }
  return(sum(log(y_tmp)))
}

sample_z_tilde_gp_dir_K_slice_sampler <- function(
  prev_z_tilde_list, y, rho_list, p_list, gp_chol, mu_list, K) {
  
  nu <- matrix(rnorm(length(y) * ncol(prev_z_tilde_list)), ncol = ncol(prev_z_tilde_list))
  nu <- gp_chol %*% nu
  
  y_tmp = rep(0, length(y))
  log_y = calc_svm_c_ll(y, p_list, prev_z_tilde_list, rho_list, mu_list, K) + 
    log(runif(1))
  theta <- runif(1, 0, 2 * pi)
  theta_b <- c(theta - 2 * pi, theta)
  
  repeat {
    next_z_tilde_list <- prev_z_tilde_list * cos(theta) + nu * sin(theta)
    next_ll <- calc_svm_c_ll(y, p_list, next_z_tilde_list, rho_list, mu_list, K)
    if (next_ll > log_y) {
      break
    }
    if (theta < 0) {
      theta_b[1] = theta
    } else {
      theta_b[2] = theta
    }
    theta <- runif(1, theta_b[1], theta_b[2])
  }
  return(next_z_tilde_list)
  
}

sample_rho_tilde_i <- function(rho_tilde_i, y_i, m_i, rho_tilde_mu) {
  next_rho_tilde_i = rho_tilde_i + rnorm(1, 0, 0.1)
  next_rho_i = exp(next_rho_tilde_i)
  log_accept_prob = 
    calc_von_mises_prob(y_i, m_i, next_rho_i, log = T) +
      dnorm(next_rho_tilde_i, rho_tilde_mu, 0.05, log = T) - 
    (calc_von_mises_prob(y_i, m_i, exp(rho_tilde_i), log = T) +
       dnorm(rho_tilde_i, rho_tilde_mu, 0.05, log = T))
  if (log(runif(1)) < log_accept_prob) {
    return(next_rho_tilde_i)
  } 
  return(rho_tilde_i)
}

sample_rho_tilde_mu <- function(rho_tilde_list) {
  post_var <- 1/5^2 + length(rho_tilde_list) / 0.05^2
  post_mean <- 1 / post_var * (sum(rho_tilde_list) / 0.05^2)
  return(rnorm(1,  post_mean, 1 / sqrt(post_var)))
}

sample_p_list <- function(prev_p_list, y, m_list, rho_list) {
  prev_log_post <- calc_svm_c_ll(y, prev_p_list, NA, rho_list, NA, m_list)
  prev_p_list <- logit(prev_p_list)
  for (i in sample(length(prev_p_list))) {
    prop_info <- PropStep(prev_p_list, i, rep(0.1, length(prev_p_list)))
    next_p_list = prop_info[1:length(prev_p_list)]
    next_log_post <- calc_svm_c_ll(y, invlogit(next_p_list), NA, rho_list, NA, m_list)
    log_accept_prob = attr(prop_info, "dbt") + 
      next_log_post - prev_log_post
    if (log(runif(1)) < log_accept_prob) {
      prev_p_list = next_p_list
      prev_log_post = next_log_post
    }
  }
  return(invlogit(prev_p_list))
}
 
generate_sim_post_pred_prob <- function(m_tilde_list, pi_param, full_cov, 
                                        mean_m_list, pred_dir, K, p_list,
                                        obs_ind_v = 1:500, pred_ind_v = 501:550) {
  cov_chol <- t(chol(full_cov[obs_ind_v, obs_ind_v]))
  pred_mu <- solve(t(cov_chol), m_tilde_list)
  pred_mu <- full_cov[pred_ind_v, obs_ind_v] %*% pred_mu
  pred_mu <- t(apply(pred_mu, 1, function(row) {row + mean_m_list}))

  cond_cov <- full_cov[pred_ind_v, pred_ind_v] - 
    full_cov[pred_ind_v, obs_ind_v] %*% solve(full_cov[obs_ind_v, obs_ind_v], full_cov[obs_ind_v, pred_ind_v])
  cond_cov_chol <- t(chol(cond_cov))

  lp_predict <- rep(0, 100)
  for (i in 1:length(lp_predict)) {
    rand_pred_z_draw <- apply(pred_mu, 2, function(col) {
      tmp <- rnorm(length(col))
      cond_cov_chol %*% tmp + col})
    tmp <- rep(0, length(pred_dir))
    for (k in 1:K) {
      pred_m_draw <- 
        mapply(function(z1, z2) {angle(z1, z2)}, 
              rand_pred_z_draw[, (k - 1) * 2 + 1],
              rand_pred_z_draw[, (k - 1) * 2 + 2])
      pi_param_list <- exp(rnorm(length(pred_dir), pi_param[k], 0.05))
      tmp = tmp + p_list[k] * calc_von_mises_prob(pred_dir, pred_m_draw, pi_param_list)
    }
    lp_predict[i] = sum(log(tmp))
  }
  return(lp_predict)
}

sample_svm <- function(y, x, sigma, l, mean_v, num_iter = 2000,
                       z_tilde_start = NULL, z_tilde_cov = F,
                       rho_tilde_i_start = NULL,
                       rho_tilde_mu_start = NULL, 
                       calc_pred_prob = F, pred_x = NULL,
                       pred_y = NULL, gp_dir_model = NULL) {
  
  cov <- as.matrix(dist(x))^2
  cov <- sigma^2 * exp(-cov / (2 * l^2)) + diag(rep(1e-9, nrow(cov)))
  cov_chol <- t(chol(cov))
  

  z_tilde_draws <- matrix(0, num_iter + 1, ncol = 2 * length(y))
  m_draws <- matrix(angle(mean_v[1], mean_v[2]), nrow = num_iter + 1, ncol = length(y))
  if (!is.null(z_tilde_start)) {
    if (!z_tilde_cov) {
      z_tilde_start = cov_chol %*% z_tilde_start
    }
    z_tilde_draws[1,] <- as.vector(z_tilde_start)
    m_draws[1,] <- apply(z_tilde_start, 1, function(start) {
      angle(start[1] + mean_v[1], start[2] + mean_v[2])
    })
  }
  
  rho_tilde_i_draws <- matrix(0, nrow = num_iter + 1, ncol = length(y))
  if (!is.null(rho_tilde_i_start)) {
    rho_tilde_i_draws[1,] = rho_tilde_i_start
  } 

  rho_tilde_mu_draws <- rep(0, num_iter + 1)
  if (!is.null(rho_tilde_mu_start)) {
    rho_tilde_mu_draws[1] = rho_tilde_mu_start
  } 
  if (calc_pred_prob) {
    full_cov <- as.matrix(dist(rbind(x, pred_x)))^2
    full_cov <- sigma^2 * exp(-full_cov / (2 * l^2)) + diag(rep(1e-9, nrow(full_cov))) 
    post_pred_prob <- matrix(0, nrow = num_iter, ncol = 100)
  }
  for (i in 1:num_iter) {
    if (i %% 100 == 1) {
      print(i)
    }
    z_tilde_draws[i + 1,] <- as.vector(
      sample_z_tilde_gp_dir_slice_sampler(
        matrix(z_tilde_draws[i,], ncol = 2), y, exp(rho_tilde_i_draws[i,]), cov_chol, mean_v))
    m_draws[i + 1,] <- sapply(1:length(y), function(j) {
      angle(z_tilde_draws[i + 1, j] + mean_v[1], z_tilde_draws[i + 1, j + length(y)] + mean_v[2])
    })
    
    data_list <- list("N" = length(y), "observed_dir" = y,
                      "trans_m_list" = as.vector(m_draws[i + 1,]))
    init_list <- list(list("pi_param_tilde_list" = (as.vector(rho_tilde_i_draws[i,]) - rho_tilde_mu_draws[i]) / 0.05,
                      "pi_param_mu" = rho_tilde_mu_draws[i]))
    gp_dir_results <- sampling(gp_dir_model, data = data_list, init = init_list,
                               chains = 1, iter = 5, warmup = 4,
                               open_progress = F, show_messages = F,
                               control = list("adapt_delta" = 0.99))
    gp_dir_info <- as.data.frame(gp_dir_results)
    rho_tilde_i_draws[i + 1,] <- 0.05 * unname(unlist(gp_dir_info[nrow(gp_dir_info), grep("pi_param_tilde_list", colnames(gp_dir_info))])) + 
      gp_dir_info[nrow(gp_dir_info), "pi_param_mu"]
    rho_tilde_mu_draws[i + 1] <- gp_dir_info[nrow(gp_dir_info), "pi_param_mu"]
    
    if (calc_pred_prob) {
      post_pred_prob[i,] <-
        generate_sim_post_pred_prob(solve(cov_chol, 
                                      matrix(z_tilde_draws[i + 1,], ncol = 2)), 
                                    rho_tilde_mu_draws[i + 1], full_cov, 
                                    mean_v, pred_y, 1, 1,
                                    obs_ind_v = 1:length(y), 
                                    pred_ind_v = length(y) + 1:length(pred_y))
    }
  }
  if (calc_pred_prob) {
    return(list(m_draws[-1,], z_tilde_draws[-1,], rho_tilde_i_draws[-1,], rho_tilde_mu_draws[-1],
                post_pred_prob))
  }
  return(list(m_draws[-1,], z_tilde_draws[-1,], rho_tilde_i_draws[-1,], rho_tilde_mu_draws[-1]))
}

