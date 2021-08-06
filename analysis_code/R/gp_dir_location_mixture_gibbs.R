library(rstan)
library(Rcpp)
library(LaplacesDemon)
library(SALTSampler)

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


create_trans_m_matrix <- function(z_tilde_list, mu_list, K) {
  trans_m_list <- matrix(0, nrow = nrow(z_tilde_list), ncol = K)
  for (k in 1:K) {
    trans_m_list[,k] <- mapply(function(tilde_z1, tilde_z2, mu_1, mu_2) {
        angle(tilde_z1 + mu_1, tilde_z2 + mu_2)
      }, 
      z_tilde_list[, 2 * (k - 1) + 1], z_tilde_list[, 2 * (k - 1) + 2],
      mu_list[2 * (k - 1) + 1], mu_list[2 * (k - 1) + 2])
  }
  return(trans_m_list)
}

calc_svm_c_ll <- function(y, p_list, z_tilde_list, rho_list, mu_list, K, comp_m_list = NULL) {
  y_tmp = rep(0, length(y))
  if (is.null(comp_m_list)) {
    comp_m_list <- create_trans_m_matrix(z_tilde_list, mu_list, K)
  }
  for (k in 1:K) {
    y_tmp = y_tmp + p_list[k] * calc_von_mises_prob(y, comp_m_list[, k], rho_list[, k])
  }
  return(sum(log(y_tmp)))
}

sample_z_tilde_gp_dir_slice_sampler <- function(
  prev_z_tilde_list, y, rho_list, gp_chol, mu_list, sample_inds = 1:length(y)) {
  
  nu <- matrix(rnorm(length(y) * ncol(prev_z_tilde_list)), ncol = 2)
  nu <- gp_chol %*% nu
  prev_m_list <- apply(prev_z_tilde_list, 1, function(row) {
    angle(row[1] + mu_list[1], row[2] + mu_list[2])
  })
  log_y <- sum(
    calc_von_mises_prob(y[sample_inds], 
                        prev_m_list[sample_inds], 
                        rho_list[sample_inds], log = T)) + log(runif(1))
  theta <- runif(1, 0, 2 * pi)
  theta_b <- c(theta - 2 * pi, theta)
  
  repeat {
    next_z_tilde_list <- prev_z_tilde_list * cos(theta) + nu * sin(theta)
    next_m_list <- apply(next_z_tilde_list, 1, function(row) {
      angle(row[1] + mu_list[1], row[2] + mu_list[2])
    })
    next_ll <- sum(
      calc_von_mises_prob(y[sample_inds], 
                          next_m_list[sample_inds], 
                          rho_list[sample_inds], log = T))
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

sample_label_list_draws <- function(p_list, y, m_list, rho_list) {
  label_probs <- sapply(1:length(p_list), function(k) {
    p_list[k] * calc_von_mises_prob(y, m_list[, k], rho_list[, k]) 
  })
  label_probs <- label_probs / rowSums(label_probs)
  apply(label_probs, 1, function(row) {
    sample(length(p_list), 1, prob = row)
  })
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

sample_rho_tilde_i_mu_svm_c_stan <- function(data_list, init_list, gp_dir_model) {
  gp_dir_results <- sampling(gp_dir_model, data = data_list, init = init_list,
                             chains = 1, iter = 10, open_progress = F, show_messages = F,
                             warmup = 9, control = list("adapt_delta" = 0.99))
  gp_dir_info <- as.data.frame(gp_dir_results)
  rho_tilde_mu_draws <- 
    data.matrix(gp_dir_info[nrow(gp_dir_info), grep("pi_param_mu", colnames(gp_dir_info))])
  rho_tilde_i_draws <- 
    0.05 * unname(unlist(gp_dir_info[nrow(gp_dir_info), grep("pi_param_tilde_list", colnames(gp_dir_info))])) + 
    rep(rho_tilde_mu_draws, each = length(data_list$observed_dir))
  return(list(rho_tilde_mu_draws, rho_tilde_i_draws))
}

sample_rho_tilde_i_mu_svm_c_stan_with_labels <- 
  function(data_list, init_list, gp_dir_model) {

  gp_dir_results <- sampling(gp_dir_model, data = data_list, init = init_list,
                             chains = 1, iter = 10, open_progress = F, show_messages = F,
                             warmup = 9, control = list("adapt_delta" = 0.99))
  gp_dir_info <- as.data.frame(gp_dir_results)
  rho_tilde_mu_draws <- 
    data.matrix(gp_dir_info[nrow(gp_dir_info), grep("pi_param_mu", colnames(gp_dir_info))])
  rho_tilde_i_draws <- 
    0.05 * unname(unlist(gp_dir_info[nrow(gp_dir_info), grep("pi_param_tilde_list", colnames(gp_dir_info))])) + 
    rep(rho_tilde_mu_draws, each = length(data_list$observed_dir))
  return(list(rho_tilde_mu_draws, rho_tilde_i_draws))
}

sample_p_list <- function(prev_p_list, y, m_list, rho_list) {
  prev_log_post <- calc_svm_c_ll(y, prev_p_list, NA, rho_list, NA, length(prev_p_list), m_list)
  prev_p_list <- logit(prev_p_list)
  for (i in sample(length(prev_p_list))) {
    prop_info <- PropStep(prev_p_list, i, rep(0.1, length(prev_p_list)))
    next_p_list = prop_info[1:length(prev_p_list)]
    next_log_post <- calc_svm_c_ll(y, invlogit(next_p_list), NA, rho_list, NA, length(prev_p_list), 
                                   m_list)
    log_accept_prob = attr(prop_info, "dbt") + 
      next_log_post - prev_log_post
    if (log(runif(1)) < log_accept_prob) {
      prev_p_list = next_p_list
      prev_log_post = next_log_post
    }
  }
  return(invlogit(prev_p_list))
}

sample_svm_c <- function(y, x, sigma, l, mean_v, K, num_iter = 2000,
                         z_tilde_start = NULL, z_tilde_cov = F,
                         rho_tilde_i_start = NULL,
                         rho_tilde_mu_start = NULL, 
                         p_list_start = NULL,
                         calc_pred_prob = F, pred_x = NULL,
                         pred_y = NULL, gp_dir_model = NULL) {
  
  cov <- as.matrix(dist(x))^2
  cov <- sigma^2 * exp(-cov / (2 * l^2)) + diag(rep(1e-9, nrow(cov)))
  cov_chol <- t(chol(cov))
  
  p_list_draws <- matrix(1 / K, num_iter + 1, ncol = K)
  if (!is.null(p_list_draws)) {
    p_list_draws[1,] <- p_list_start
  }
  
  z_tilde_draws <- matrix(0, num_iter + 1, ncol = 2 * K * length(y))
  start_angles <- rep(0, K)
  for (k in 1:K) {
    start_angles[k] = angle(mean_v[2 * (k - 1) + 1], mean_v[2 * (k - 1) + 2])
  }
  m_draws <- matrix(rep(start_angles, each = length(y)), 
                    nrow = num_iter + 1, ncol = K * length(y), byrow = T)
  if (!is.null(z_tilde_start)) {
    if (!z_tilde_cov) {
      z_tilde_start = cov_chol %*% z_tilde_start
    }
    z_tilde_draws[1,] <- as.vector(z_tilde_start)
    m_draws[1,] <- 
        as.vector(create_trans_m_matrix(z_tilde_start, mean_v, K))
  }
  
  rho_tilde_i_draws <- matrix(0, nrow = num_iter + 1, ncol = K * length(y))
  if (!is.null(rho_tilde_i_start)) {
    rho_tilde_i_draws[1,] = as.vector(rho_tilde_i_start)
  } 
  
  rho_tilde_mu_draws <- matrix(0, nrow = num_iter + 1, ncol = K)
  if (!is.null(rho_tilde_mu_start)) {
    rho_tilde_mu_draws[1,] = rho_tilde_mu_start
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
    z_tilde_tmp <- sample_z_tilde_gp_dir_K_slice_sampler(
      matrix(z_tilde_draws[i,], ncol = 2 * K), y, matrix(exp(rho_tilde_i_draws[i,]), ncol = K), 
      p_list_draws[i,], cov_chol, mean_v, K)
    z_tilde_draws[i + 1,] <- as.vector(z_tilde_tmp)
    m_draws[i + 1,] <-  as.vector(create_trans_m_matrix(z_tilde_tmp, mean_v, K))

    data_list <- list("N" = length(y), "K" = K, "observed_dir" = y, "p_list" = p_list_draws[i,],
                      "trans_m_list" = matrix(m_draws[i + 1,], ncol = K))
    init_list <- list(list("pi_param_tilde_list" = 
                             matrix((rho_tilde_i_draws[i,] - rep(rho_tilde_mu_draws[i,], each = length(y))) / 0.05, ncol = K),
                           "pi_param_mu" = rho_tilde_mu_draws[i,]))
    gp_dir_results <- sampling(gp_dir_model, data = data_list, init = init_list,
                               chains = 1, iter = 10, open_progress = F, show_messages = F,
                               warmup = 9, control = list("adapt_delta" = 0.99))
    gp_dir_info <- as.data.frame(gp_dir_results)
    rho_tilde_mu_draws[i + 1,] <- 
      data.matrix(gp_dir_info[nrow(gp_dir_info), grep("pi_param_mu", colnames(gp_dir_info))])
    rho_tilde_i_draws[i + 1,] <- 
      0.05 * unname(unlist(gp_dir_info[nrow(gp_dir_info), grep("pi_param_tilde_list", colnames(gp_dir_info))])) + 
      rep(rho_tilde_mu_draws[i + 1,], each = length(y))
    
    p_list_draws[i + 1,] <- sample_p_list(p_list_draws[i,], y, matrix(m_draws[i + 1,], ncol = K),
                                          matrix(exp(rho_tilde_i_draws[i + 1,]), ncol = K))

    if (calc_pred_prob) {
      post_pred_prob[i,] <-
        generate_sim_post_pred_prob(solve(cov_chol, matrix(z_tilde_draws[i + 1,], ncol = 2 * K)), 
                                    rho_tilde_mu_draws[i + 1,], full_cov, 
                                    mean_v, pred_y, K, p_list_draws[i + 1,],
                                    1:length(y), length(y) + 1:length(pred_y))
    }
  }
  if (calc_pred_prob) {
    return(list(p_list_draws[-1,], m_draws[-1,], z_tilde_draws[-1,], rho_tilde_i_draws[-1,], rho_tilde_mu_draws[-1,],
                post_pred_prob))
  }
  return(list(p_list_draws[-1,], m_draws[-1,], z_tilde_draws[-1,], rho_tilde_i_draws[-1,], rho_tilde_mu_draws[-1,]))
}

sample_svm_c_labels <- function(y, x, sigma, l, mean_v, K, num_iter = 2000,
                         z_tilde_start = NULL, z_tilde_cov = F,
                         rho_tilde_i_start = NULL,
                         rho_tilde_mu_start = NULL, 
                         p_list_start = NULL, label_list_start = NULL,
                         calc_pred_prob = F, pred_x = NULL,
                         pred_y = NULL, gp_dir_model = NULL) {
  
  cov <- as.matrix(dist(x))^2
  cov <- sigma^2 * exp(-cov / (2 * l^2)) + diag(rep(1e-9, nrow(cov)))
  cov_chol <- t(chol(cov))
  
  p_list_draws <- matrix(1 / K, nrow = num_iter + 1, ncol = K)
  label_list_draws <- matrix(sample(K, length(y), replace = T), 
                             nrow = num_iter + 1, ncol = length(y), 
                             byrow = T)
  if (!is.null(p_list_start)) {
    p_list_draws[1,] <- p_list_start
    label_list_draws[1,] <- sample(K, length(y), 
                                   replace = T, prob = p_list_start)
  }
  if (!is.null(label_list_start)) {
    label_list_draws[1,] <- label_list_start
  }  

  z_tilde_draws <- matrix(0, num_iter + 1, ncol = 2 * K * length(y))
  start_angles <- rep(0, K)
  for (k in 1:K) {
    start_angles[k] = angle(mean_v[2 * (k - 1) + 1], mean_v[2 * (k - 1) + 2])
  }
  m_draws <- matrix(rep(start_angles, each = length(y)), 
                    nrow = num_iter + 1, ncol = K * length(y), byrow = T)
  if (!is.null(z_tilde_start)) {
    if (!z_tilde_cov) {
      z_tilde_start = cov_chol %*% z_tilde_start
    }
    z_tilde_draws[1,] <- as.vector(z_tilde_start)
    m_draws[1,] <- 
        as.vector(create_trans_m_matrix(z_tilde_start, mean_v, K))
  }
  
  rho_tilde_i_draws <- matrix(0, nrow = num_iter + 1, ncol = K * length(y))
  if (!is.null(rho_tilde_i_start)) {
    rho_tilde_i_draws[1,] = as.vector(rho_tilde_i_start)
  } 
  
  rho_tilde_mu_draws <- matrix(0, nrow = num_iter + 1, ncol = K)
  if (!is.null(rho_tilde_mu_start)) {
    rho_tilde_mu_draws[1,] = rho_tilde_mu_start
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
    
    #Update labels
    label_list_draws[i + 1,] <- sample_label_list_draws(p_list_draws[i,], y,  
                                                        matrix(m_draws[i,], ncol = K),
                                                        matrix(exp(rho_tilde_i_draws[i,]), ncol = K))
    
    label_counts <- rep(0, K)
    label_count_tmp <- table(label_list_draws[i + 1,])
    label_counts[as.numeric(names(label_count_tmp))] <-
      label_count_tmp
    p_list_draws[i + 1,] <- rdirichlet(1, label_counts + rep(1, K))
    
    #Update m_draws and r_draws
    z_tilde_tmp <- do.call(cbind, lapply(1:K, function(k) {
      cluster_ind <- which(label_list_draws[i + 1,] == k)
      sample_z_tilde_gp_dir_slice_sampler(
        matrix(z_tilde_draws[i, 2 * (k - 1) * length(y) + 1:(2 * length(y))], ncol = 2), 
        y,
        exp(rho_tilde_i_draws[i, (k - 1) * length(y) + 1:length(y)]), 
        cov_chol, mean_v[2 * (k - 1) + 1:2],
        cluster_ind)      
    }))
    z_tilde_draws[i + 1,] <- as.vector(z_tilde_tmp)
    m_draws[i + 1,] <-  as.vector(create_trans_m_matrix(z_tilde_tmp, mean_v, K))

    #Update concentration parameters
    data_list <- list("N" = length(y), "K" = K, "observed_dir" = y, "p_list" = p_list_draws[i,],
                      "trans_m_list" = matrix(m_draws[i + 1,], ncol = K),
                      "labels" = label_list_draws[i + 1,])
    init_list <- list(list("pi_param_tilde_list" =
                             matrix((rho_tilde_i_draws[i,] - rep(rho_tilde_mu_draws[i,], each = length(y))) / 0.05, ncol = K),
                           "pi_param_mu" = rho_tilde_mu_draws[i,]))
    rho_tilde_i_mu_info <- sample_rho_tilde_i_mu_svm_c_stan_with_labels(data_list, init_list, gp_dir_model)
    rho_tilde_mu_draws[i + 1,] <- rho_tilde_i_mu_info[[1]]
    rho_tilde_i_draws[i + 1,] <- rho_tilde_i_mu_info[[2]]

    if (calc_pred_prob) {
      post_pred_prob[i,] <-
        generate_sim_post_pred_prob(solve(cov_chol, matrix(z_tilde_draws[i + 1,], ncol = 2 * K)), 
                                    rho_tilde_mu_draws[i + 1,], full_cov, 
                                    mean_v, pred_y, K, p_list_draws[i + 1,], 
                                    1:length(y), length(y) + 1:length(pred_y))
    }
  }
  if (calc_pred_prob) {
    return(list(p_list_draws[-1,], m_draws[-1,], z_tilde_draws[-1,], rho_tilde_i_draws[-1,], rho_tilde_mu_draws[-1,],
                post_pred_prob, label_list_draws[-1,]))
  }
  return(list(p_list_draws[-1,], m_draws[-1,], z_tilde_draws[-1,], rho_tilde_i_draws[-1,], rho_tilde_mu_draws[-1,], label_list_draws[-1,]))
}


