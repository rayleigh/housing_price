library(matrixStats)
library(MCMCpack)

#Higher dimensional transformation
calc_von_mises_fisher_prob <- function(y, mean_m, rho, log = T) {
  circ_dim = ncol(y)
  log_prob = (circ_dim / 2 - 1) * log(rho)
  log_prob = log_prob - log(besselI(rho, circ_dim / 2 - 1)) - 
    circ_dim / 2 * log(2 * pi)
  # log_prob = log_prob * nrow(y)
  log_prob = log_prob + rho * rowSums(y * mean_m)
  if (log == F) {
    return(exp(log_prob))
  }
  return(log_prob)
}

calc_expected_log_likelihood_k <- function(
  y, p, m_k_m, tilde_p_k_m, label_v, k, 
  tilde_z_k_m = NULL, mu_tilde_p_k_val = NULL, 
  centered = T, noncentered_tilde_p_k_m = NULL) {
  
  p_k_m <- exp(tilde_p_k_m)
  log_ll = 0
  if (!is.null(mu_tilde_p_k_val)) {
    log_ll = log_ll + 
      sum(dnorm(mu_tilde_p_k_val, 0, 5, log = T))
    if (centered) {
      log_ll = log_ll + 
        sum(dnorm(tilde_p_k_m, mu_tilde_p_k_val, 0.05, log = T))
    } else {
      log_ll = log_ll + 
        sum(dnorm(noncentered_tilde_p_k_m, 0, 1, log = T))
    }
      
  }
  if (!is.null(tilde_z_k_m)) {
    log_ll = log_ll + sum(dnorm(tilde_z_k_m, log = T))
  }
  # log_ll = sum(dnorm(mu_tilde_p_k, 0, 5, log = T)) +
  #   sum(dnorm(tilde_z_k_m, log = T))
  
  log_ll = log_ll + sum(
    (label_v == k) * (
      log(p[k]) +
        calc_von_mises_fisher_prob(y, m_k_m, p_k_m, log = T))) 

  return(log_ll)
}

grad_update_concentration_hmc_with_label <- 
  function(next_tilde_p_k, next_p_k, m_k_m, y_m, mu_tilde_p_k_val, label_v, k) {
  
    # ll_grad <- ll_grad - (
    #   besselI(next_p_k, d / 2 - 1) * (d / 2 - 1) / (next_p_k) + 
    #     besselI(next_p_k, d / 2)) / 
    #   besselI(next_p_k, d / 2 - 1)
    
    # grad_tilde_p_k <- r_k_m[, k] * next_p_k * ll_grad
    
  d = ncol(y_m)
  ll_grad <- rowSums(m_k_m * y_m) - 
    besselI(next_p_k, d / 2) / besselI(next_p_k, d / 2 - 1)
  grad_tilde_p_k <- next_p_k * ll_grad
  grad_tilde_p_k <- grad_tilde_p_k * (label_v == k)
  # grad_tilde_p_k = grad_tilde_p_k - r_k_m[, k] / (0.1^2) * (next_tilde_p_k - mu_tilde_p_k[k])
  grad_tilde_p_k = grad_tilde_p_k - 1 / (0.05^2) * (next_tilde_p_k - mu_tilde_p_k_val)
  return(grad_tilde_p_k)
}

grad_update_hierarchical_concentration_hmc_with_label <-
  function(next_tilde_p_k, mu_tilde_p_k_val) {
    
  return(-1 / (0.05^2) * sum(next_tilde_p_k - mu_tilde_p_k_val) -
           mu_tilde_p_k_val / 5^2) 
  }


grad_update_concentration_hmc_for_k <- 
  function(tilde_p_k_m, mu_tilde_p_k_val, p, m_k_m, y_m, label_v, k,
           hmc_epsilon, hmc_steps, hmc_sd) {
    
    momentum_v <- rnorm(length(tilde_p_k_m) + 1, 0, sd = hmc_sd)
    next_tilde_p_k_m <- tilde_p_k_m
    next_mu_tilde_p_k_val <- mu_tilde_p_k_val
    
    next_momentum_v <- momentum_v + hmc_epsilon / 2 * c(
      grad_update_concentration_hmc_with_label(
        next_tilde_p_k_m, exp(next_tilde_p_k_m), m_k_m, y_m, 
        next_mu_tilde_p_k_val, label_v, k),
      grad_update_hierarchical_concentration_hmc_with_label(
        next_tilde_p_k_m, next_mu_tilde_p_k_val))
    next_tilde_p_k_m <- next_tilde_p_k_m + 
      hmc_epsilon / hmc_sd^2 * next_momentum_v[2:length(next_momentum_v) - 1]
    next_mu_tilde_p_k_val <- next_mu_tilde_p_k_val + 
      hmc_epsilon / hmc_sd^2 * next_momentum_v[length(next_momentum_v)]
    if (hmc_steps >= 2) {
      for (i in 2:hmc_steps) {
        next_momentum_v <- next_momentum_v + hmc_epsilon * c(
          grad_update_concentration_hmc_with_label(
            next_tilde_p_k_m, exp(next_tilde_p_k_m), m_k_m, y_m, 
            next_mu_tilde_p_k_val, label_v, k),
          grad_update_hierarchical_concentration_hmc_with_label(
            next_tilde_p_k_m, next_mu_tilde_p_k_val))
        next_tilde_p_k_m <- next_tilde_p_k_m + 
          hmc_epsilon / hmc_sd^2 * next_momentum_v[2:length(next_momentum_v) - 1]
        next_mu_tilde_p_k_val <- next_mu_tilde_p_k_val + 
          hmc_epsilon / hmc_sd^2 * next_momentum_v[length(next_momentum_v)]
        # print(next_momentum_v)
        # print(next_tilde_p_k_m)
        # print(next_mu_tilde_p_k_val)
      }  
    }
    next_momentum_v <- next_momentum_v + hmc_epsilon / 2 * c(
      grad_update_concentration_hmc_with_label(
        next_tilde_p_k_m, exp(next_tilde_p_k_m), m_k_m, y_m, 
        next_mu_tilde_p_k_val, label_v, k),
      grad_update_hierarchical_concentration_hmc_with_label(
        next_tilde_p_k_m, next_mu_tilde_p_k_val))
    
    momentum_prob <- 
      sum(momentum_v^2) / (2 * hmc_sd^2) -
      sum(next_momentum_v^2) / (2 * hmc_sd^2)
    # y, p, m_k_m, tilde_p_k_m, label_v, k, 
    # tilde_z_k_m = NULL, mu_tilde_p_k_val = NULL
    U <- calc_expected_log_likelihood_k(
      y_m, p, m_k_m, next_tilde_p_k_m, label_v, k, 
      mu_tilde_p_k_val = next_mu_tilde_p_k_val) -
      calc_expected_log_likelihood_k(
        y_m, p, m_k_m, tilde_p_k_m, label_v, k,
        mu_tilde_p_k_val = mu_tilde_p_k_val)
    if (log(runif(1)) < U + momentum_prob) {
      return(c(next_tilde_p_k_m, next_mu_tilde_p_k_val))
    }
    return(c(tilde_p_k_m, mu_tilde_p_k_val))
  }

ll_grad_with_label <- 
  function(next_p_k, m_k_m, y_m, label_v, k) {
    
    d = ncol(y_m)
    ll_grad <- rowSums(m_k_m * y_m) - 
      besselI(next_p_k, d / 2) / besselI(next_p_k, d / 2 - 1)
    grad_tilde_p_k <- next_p_k * ll_grad
    grad_tilde_p_k <- grad_tilde_p_k * (label_v == k)
    return(grad_tilde_p_k)
}

grad_update_concentration_hmc_noncentered_with_label <- 
  function(
    next_noncentered_tilde_p_k, 
    next_p_k, m_k_m, y_m, label_v, k) {
    
    # ll_grad <- ll_grad - (
    #   besselI(next_p_k, d / 2 - 1) * (d / 2 - 1) / (next_p_k) + 
    #     besselI(next_p_k, d / 2)) / 
    #   besselI(next_p_k, d / 2 - 1)
    
    # grad_tilde_p_k <- r_k_m[, k] * next_p_k * ll_grad
    
    grad_tilde_p_k <- ll_grad_with_label(next_p_k, m_k_m, y_m, label_v, k)
    grad_tilde_p_k <- 0.05 * grad_tilde_p_k - next_noncentered_tilde_p_k
    # grad_tilde_p_k = grad_tilde_p_k - r_k_m[, k] / (0.1^2) * (next_tilde_p_k - mu_tilde_p_k[k])
    # grad_tilde_p_k = grad_tilde_p_k - 1 / (0.05^2) * (next_tilde_p_k - mu_tilde_p_k_val)
    return(grad_tilde_p_k)
  }

grad_update_hierarchical_concentration_noncentered_hmc_with_label <-
  function(next_p_k, m_k_m, y_m, mu_tilde_p_k_val, label_v, k) {
    
   return(sum(ll_grad_with_label(next_p_k, m_k_m, y_m, label_v, k)) - 
            mu_tilde_p_k_val / 5^2)
}

grad_update_noncentered_concentration_hmc_for_k <- 
  function(tilde_p_k_m, mu_tilde_p_k_val, p, m_k_m, y_m, label_v, k,
           hmc_epsilon, hmc_steps, hmc_sd) {
    
    momentum_v <- rnorm(length(tilde_p_k_m) + 1, 0, sd = hmc_sd)
    next_tilde_p_k_m <- tilde_p_k_m
    next_mu_tilde_p_k_val <- mu_tilde_p_k_val
    next_noncentered_tilde_p_k_m <- 
      (tilde_p_k_m - mu_tilde_p_k_val) / 0.05
    
    next_momentum_v <- momentum_v + hmc_epsilon / 2 * c(
      grad_update_concentration_hmc_noncentered_with_label(
        next_noncentered_tilde_p_k_m, exp(next_tilde_p_k_m), m_k_m, y_m, 
        label_v, k),
      grad_update_hierarchical_concentration_noncentered_hmc_with_label(
        exp(next_tilde_p_k_m), m_k_m, y_m, next_mu_tilde_p_k_val,
        label_v, k))
    next_noncentered_tilde_p_k_m <- next_noncentered_tilde_p_k_m + 
      hmc_epsilon / hmc_sd^2 * next_momentum_v[2:length(next_momentum_v) - 1]
    next_mu_tilde_p_k_val <- next_mu_tilde_p_k_val + 
      hmc_epsilon / hmc_sd^2 * next_momentum_v[length(next_momentum_v)]
    next_tilde_p_k_m <- 0.05 * next_noncentered_tilde_p_k_m + 
      next_mu_tilde_p_k_val
    if (hmc_steps >= 2) {
      for (i in 2:hmc_steps) {
        next_momentum_v <- next_momentum_v + hmc_epsilon * c(
          grad_update_concentration_hmc_noncentered_with_label(
            next_noncentered_tilde_p_k_m, exp(next_tilde_p_k_m), m_k_m, y_m, 
            label_v, k),
          grad_update_hierarchical_concentration_noncentered_hmc_with_label(
            exp(next_tilde_p_k_m), m_k_m, y_m, next_mu_tilde_p_k_val,
            label_v, k))
        next_noncentered_tilde_p_k_m <- next_noncentered_tilde_p_k_m + 
          hmc_epsilon / hmc_sd^2 * next_momentum_v[2:length(next_momentum_v) - 1]
        next_mu_tilde_p_k_val <- next_mu_tilde_p_k_val + 
          hmc_epsilon / hmc_sd^2 * next_momentum_v[length(next_momentum_v)]
        next_tilde_p_k_m <- 0.05 * next_noncentered_tilde_p_k_m + 
          next_mu_tilde_p_k_val
        # print(next_momentum_v)
        # print(next_tilde_p_k_m)
        # print(next_mu_tilde_p_k_val)
      }  
    }
    next_momentum_v <- next_momentum_v + hmc_epsilon / 2 * c(
      grad_update_concentration_hmc_noncentered_with_label(
        next_noncentered_tilde_p_k_m, exp(next_tilde_p_k_m), m_k_m, y_m, 
        label_v, k),
      grad_update_hierarchical_concentration_noncentered_hmc_with_label(
        exp(next_tilde_p_k_m), m_k_m, y_m, next_mu_tilde_p_k_val,
        label_v, k))
    
    momentum_prob <- 
      sum(momentum_v^2) / (2 * hmc_sd^2) -
      sum(next_momentum_v^2) / (2 * hmc_sd^2)
    # y, p, m_k_m, tilde_p_k_m, label_v, k, 
    # tilde_z_k_m = NULL, mu_tilde_p_k_val = NULL
    U <- calc_expected_log_likelihood_k(
      y_m, p, m_k_m, next_tilde_p_k_m, label_v, k, 
      mu_tilde_p_k_val = next_mu_tilde_p_k_val, 
      centered = F, noncentered_tilde_p_k_m = 
        next_noncentered_tilde_p_k_m) -
      calc_expected_log_likelihood_k(
        y_m, p, m_k_m, tilde_p_k_m, label_v, k,
        mu_tilde_p_k_val = mu_tilde_p_k_val,
        centered = F, noncentered_tilde_p_k_m = 
          (tilde_p_k_m - mu_tilde_p_k_val) / 0.05)
    if (log(runif(1)) < U + momentum_prob) {
      return(c(next_tilde_p_k_m, next_mu_tilde_p_k_val))
    }
    return(c(tilde_p_k_m, mu_tilde_p_k_val))
  }


ess_sample_z_k_m_dim <- function(
  y, p, m_k_m, z_k_m, d,
  tilde_p_k_m, label_v, k, cov_chol, mu_z_k) {
  
  v <- cov_chol %*% matrix(rnorm(nrow(z_k_m)),
                           nrow = nrow(z_k_m),
                           1)
  cutoff_log_ll = calc_expected_log_likelihood_k(
    y, p, m_k_m, tilde_p_k_m, label_v, k) + log(runif(1))
  angle = runif(1, 0, 2 * pi)
  min_angle = angle - 2 * pi
  max_angle = angle
  next_z_k_m <- z_k_m
  
  repeat {
    next_z_k_m[, d] <- 
      z_k_m[, d] * cos(angle) + v * sin(angle)
    tmp <- sweep(next_z_k_m, 2, mu_z_k, "+")
    # next_z_k_m <- 
    #   sweep(z_k_m * cos(angle) + v * sin(angle), 2, 
    # mu_z_k, "+")
    next_m_k_m <- tmp / 
      apply(tmp, 1, norm, type = "2")
    next_log_ll = calc_expected_log_likelihood_k(
      y, p, next_m_k_m, tilde_p_k_m, label_v, k)
    if (next_log_ll > cutoff_log_ll) {
      break
    }
    if (angle < 0) {
      min_angle = angle
    } else {
      max_angle = angle
    }
    angle = runif(1, min_angle, max_angle)
    if (abs(angle) < .Machine$double.neg.eps) {
      next_z_k_m <- z_k_m
      break
    }
  }
  return(next_z_k_m)
}


ess_sample_z_k_m <- function(
  y, p, m_k_m, z_k_m, 
  tilde_p_k_m, label_v, k, cov_chol, mu_z_k) {
  
  v <- cov_chol %*% matrix(rnorm(length(z_k_m)),
                           nrow = nrow(z_k_m),
                           ncol = ncol(z_k_m))
  cutoff_log_ll = calc_expected_log_likelihood_k(
    y, p, m_k_m, tilde_p_k_m, label_v, k) + log(runif(1))
  angle = runif(1, 0, 2 * pi)
  min_angle = angle - 2 * pi
  max_angle = angle
  
  repeat {
    next_z_k_m <- 
      z_k_m * cos(angle) + v * sin(angle)
    tmp <- sweep(next_z_k_m, 2, mu_z_k, "+")
    # next_z_k_m <- 
    #   sweep(z_k_m * cos(angle) + v * sin(angle), 2, 
            # mu_z_k, "+")
    next_m_k_m <- tmp / 
      apply(tmp, 1, norm, type = "2")
    next_log_ll = calc_expected_log_likelihood_k(
      y, p, next_m_k_m, tilde_p_k_m, label_v, k)
    if (next_log_ll > cutoff_log_ll) {
      break
    }
    if (angle < 0) {
      min_angle = angle
    } else {
      max_angle = angle
    }
    angle = runif(1, min_angle, max_angle)
    if (abs(angle) < .Machine$double.neg.eps) {
      next_z_k_m <- z_k_m
      break
    }
  }
  return(next_z_k_m)
}

sample_label_v <- function(y, p, m_k_m_array, tilde_p_k_m, K) {
  
  p_k_m <- exp(tilde_p_k_m)
  label_prob <- sapply(1:K, function(k) {
    log(p[k]) +
      calc_von_mises_fisher_prob(y, m_k_m_array[,,k], p_k_m[,k], log = T)
  })
  return(apply(label_prob, 1, function(row) {
    sample(K, 1, prob = exp(row - logSumExp(row)))
  }))
}

vec_log_sum_exp <- function(a, b) {
  
  pmax(a, b) + log(1 + exp(pmin(a, b) - pmax(a, b)))
  
}

generate_sim_post_pred_prob_higher_dim <- function(
  z_centered_draws, pi_param, full_cov, 
  mean_m_list, pred_dir, K, p_list,
  obs_ind_v = 1:500, pred_ind_v = 501:550) {
  
  cov_chol <- t(chol(full_cov[obs_ind_v, obs_ind_v]))
  mu_chol <- t(solve(cov_chol, full_cov[obs_ind_v, pred_ind_v]))
  # pred_mu <- solve(t(cov_chol), z_centered_draws)
  # pred_mu <- full_cov[pred_ind_v, obs_ind_v] %*% pred_mu
  
  
  cond_cov <- full_cov[pred_ind_v, pred_ind_v] - 
    full_cov[pred_ind_v, obs_ind_v] %*% 
    solve(full_cov[obs_ind_v, obs_ind_v], full_cov[obs_ind_v, pred_ind_v])
  cond_cov_chol <- t(chol(cond_cov))
  
  lp_predict <- rep(0, 100)
  for (i in 1:length(lp_predict)) {
    # rand_pred_z_draw <- apply(pred_mu, 2, function(col) {
    #   tmp <- rnorm(length(col))
    #   cond_cov_chol %*% tmp + col})
    tmp <- rep(-Inf, nrow(pred_dir))
    for (k in 1:K) {
      pred_mu <- mu_chol %*% solve(cov_chol, z_centered_draws[,,k])
      pred_mu <- sweep(pred_mu, 2, mean_m_list[k,], "+")
      pred_m_draw <- matrix(rnorm(length(pred_dir)), 
                            nrow = nrow(pred_dir), ncol = ncol(pred_dir))
      pred_m_draw <- cond_cov_chol %*% pred_m_draw + pred_mu
      # pred_mu <- t(apply(pred_mu, 1, function(row) {row + mean_m_list}))
      pred_m_draw <- pred_m_draw / apply(pred_m_draw, 1, norm, type = "2")
        # mapply(function(z1, z2) {angle(z1, z2)}, 
        #        rand_pred_z_draw[, (k - 1) * 2 + 1],
        #        rand_pred_z_draw[, (k - 1) * 2 + 2])
      pi_param_list <- exp(rnorm(nrow(pred_dir), pi_param[k], 0.05))
    
      tmp = vec_log_sum_exp(tmp, log(p_list[k]) + 
        calc_von_mises_fisher_prob(pred_dir, pred_m_draw, pi_param_list))
    }
    lp_predict[i] = sum(tmp)
  }
  return(lp_predict)
}

init_data <- function(
  y, x, K, mu_z_k, model_sigma, model_l, total_iter,
  p_m_inits, z_k_m_inits, m_k_m_draws_inits,
  tilde_p_k_m_draws_inits,
  mu_tilde_p_k, mu_tilde_p_k_m_inits,
  label_inits, y_pred, x_pred,
  matern_cov, v, tilde_z = T) {
  
  if (!is.null(x_pred)) {
    model_x <- rbind(x, x_pred)
  }
  
  if (matern_cov) {
    if (v == 5/2)
      cov <- matern_kernel_five_halves(model_x, model_l)
    if (v == 3/2)
      cov <- matern_kernel_three_halves(model_x, model_l)
  } else {
    cov <- as.matrix(dist(model_x))^2
    cov <- model_sigma^2 * exp(-cov / (2 * model_l^2)) 
  }
  cov <- cov + diag(rep(1e-9, nrow(cov)))
  cov_chol <- t(chol(cov))
  
  if (!is.null(p_m_inits)) {
    p_m <- matrix(p_m_inits, ncol = K, nrow = total_iter, byrow = T)
  } else {
    p_m <- matrix(1 / K, ncol = K, nrow = total_iter)
  }
  
  if (!is.null(z_k_m_inits)) {
    z_k_m <- z_k_m_inits
    if (tilde_z) {
      for (k in 1:K) {
        # z_k_m[,,k] <- sweep(cov_chol[1:nrow(y), 1:nrow(y)] %*% 
        #                       tilde_z_k_m[,,k], 2, mu_z_k[k,], "+")
        z_k_m[,,k] <- cov_chol[1:nrow(y), 1:nrow(y)] %*%
          z_k_m[,,k]
      }
    }
    m_k_m <- z_k_m
  } else {
    tilde_z_k_m <- array(0, dim = c(nrow(y), ncol(y), K))
    z_k_m <- tilde_z_k_m
    m_k_m <- tilde_z_k_m
    for (k in 1:K) {
      # z_k_m[,,k] <- sweep(cov_chol[1:nrow(y), 1:nrow(y)] %*% 
      #                       tilde_z_k_m[,,k], 2, mu_z_k[k,], "+")
      z_k_m[,,k] <- cov_chol[1:nrow(y), 1:nrow(y)] %*%
                            tilde_z_k_m[,,k]
    }
  }
  for (k in 1:K) {
    tmp <- sweep(z_k_m[,,k], 2, mu_z_k[k,], "+")
    m_k_m[,,k] <- tmp / apply(tmp, 1, norm, type = "2")
    # m_k_m[,,k] <- z_k_m[,,k] / 
    #   apply(z_k_m[,,k], 1, norm, type = "2")
  }
  m_k_m_draws <- lapply(1:K, function(k) {
    matrix(m_k_m[,,k], nrow = total_iter, ncol = 
             ncol(y) * nrow(y))
  })
  names(m_k_m_draws) <- sapply(1:K, function(k) {
    paste("m_k_m", k, sep = "_")
  })
  
  if (!is.null(tilde_p_k_m_draws_inits)) {
    tilde_p_k_m_draws <- matrix(
      tilde_p_k_m_draws_inits, nrow = total_iter, ncol = K * nrow(y), byrow = T)
  } else {
    tilde_p_k_m_draws <- matrix(0, nrow = total_iter, ncol = K * nrow(y))
  }
  
  if (!is.null(mu_tilde_p_k_m_inits)) {
    mu_tilde_p_k_m <- matrix(mu_tilde_p_k_m_inits, nrow = total_iter, ncol = K, byrow = T)
  } else {
    mu_tilde_p_k_m <- matrix(mu_tilde_p_k, nrow = total_iter, ncol = K)
  }
  
  if (!is.null(label_inits)) {
    label_m <- matrix(label_inits, nrow = total_iter, ncol = nrow(y), byrow = T)
  } else {
    label_m <- matrix(sample(K, nrow(y), T), nrow = total_iter, ncol = nrow(y), byrow = T)
  }
  
  post_pred_prob <- matrix(0, nrow = total_iter, ncol = 100)
  return(list(p_m, z_k_m, m_k_m_draws, tilde_p_k_m_draws, 
              mu_tilde_p_k_m, label_m, post_pred_prob, cov))
}

sample_svm_c_higher_dim <- function(
  y, x, mu_z_k, model_sigma, model_l, K,
  mu_tilde_p_k = rep(0, K),
  matern_cov = F, v = 5/2,
  num_iter = 5000, start_iter = 1000, keep_iter = 1,
  hmc_epsilon = 1e-4, hmc_steps = 5, hmc_sd = 10,
  y_pred = NULL, x_pred = NULL, p_m_inits = NULL,
  z_k_m_inits = NULL, m_k_m_draws_inits = NULL,
  tilde_p_k_m_draws_inits = NULL, mu_tilde_p_k_m_inits = NULL,
  label_m_inits = NULL) {
  
  total_iter = (num_iter - start_iter) %/% keep_iter
  init_info <- init_data(
    y, x, K, mu_z_k, model_sigma, model_l, total_iter,
    p_m_inits, z_k_m_inits, m_k_m_draws_inits,
    tilde_p_k_m_draws_inits,
    mu_tilde_p_k, mu_tilde_p_k_m_inits,
    label_m_inits, y_pred, x_pred,
    matern_cov, v)

  # list(p_m, z_k_m, m_k_m_draws, tilde_p_k_m_draws, 
  #      mu_tilde_p_k_m, label_m, post_pred_prob, cov_chol)
  
  p_m <- init_info[[1]]
  p_prob <- p_m[1,]
  z_k_m <- init_info[[2]]
  m_k_m_draws <- init_info[[3]]
  m_k_m <- array(unlist(lapply(m_k_m_draws, function(draw) {
    matrix(draw[1,], ncol = ncol(y))
  })), c(nrow(y), ncol(y), K))
  tilde_p_k_m_draws <- init_info[[4]]
  tilde_p_k_m <- matrix(tilde_p_k_m_draws[1,], ncol = K)
  mu_tilde_p_k_m <- init_info[[5]]
  mu_tilde_p_k <- mu_tilde_p_k_m[1,]
  label_m <- init_info[[6]]
  label_v <- label_m[1,]
  post_pred_prob <- init_info[[7]]
  full_cov <- init_info[[8]]
  cov_chol <- t(chol(full_cov[1:nrow(x), 1:nrow(x)]))
  
  for (i in 1:num_iter) {
    if (i %% 100 == 0) {
      print(i)
    }
    
    if (K > 1) {
      label_v <- sample_label_v(y, p_prob, m_k_m, tilde_p_k_m, K)
      label_count <- rep(0, K)
      tmp <- table(label_v)
      label_count[as.numeric(names(tmp))] <- tmp
      p_prob <- rdirichlet(1, rep(1, K) + label_count)
    }
    print("Done sample")
    
    for (k in 1:K) {
      for (d in 1:ncol(z_k_m[,,k])) {
        z_k_m[,,k] <- ess_sample_z_k_m_dim(
          y, p_prob, m_k_m[,, k], z_k_m[,,k], d,
          tilde_p_k_m[,k], label_v, k, cov_chol, mu_z_k[k,])
        tmp <- sweep(z_k_m[,,k], 2, mu_z_k[k,], "+")
        m_k_m[,,k] <- tmp / apply(tmp, 1, norm, type = "2")
      }
    }
    print("Done tilde_z_k_m")
    
    for (k in 1:K) {
      p_tilde_update_info <- grad_update_noncentered_concentration_hmc_for_k(
        tilde_p_k_m[,k], mu_tilde_p_k[k], p_prob, m_k_m[,,k], y, label_v, k,
        hmc_epsilon, hmc_steps, hmc_sd)
      tilde_p_k_m[,k] <- p_tilde_update_info[2:length(p_tilde_update_info) - 1]
      mu_tilde_p_k[k] <- p_tilde_update_info[length(p_tilde_update_info)]
    }
    print("Done tilde_p_k_m")
    
    if (i > start_iter && ((i - start_iter) %% keep_iter == 0)) {
      draw_ind <- (i - start_iter) %/% keep_iter
      for (k in 1:K) {
        m_k_m_draws[[k]][draw_ind,] <- as.vector(m_k_m[,,k])
      }
      tilde_p_k_m_draws[draw_ind,] <- as.vector(tilde_p_k_m)
      mu_tilde_p_k_m[draw_ind,] <- mu_tilde_p_k
      label_m[draw_ind,] <- label_v
      p_m[draw_ind,] <- p_prob
      if (!is.null(x_pred)) {
        # z_centered_draws, pi_param, full_cov, 
        # mean_m_list, pred_dir, K, p_list,
        # obs_ind_v = 1:500, pred_ind_v = 501:550
        post_pred_prob[draw_ind,] <- generate_sim_post_pred_prob_higher_dim(
          z_k_m, mu_tilde_p_k, full_cov, 
          mu_z_k, y_pred, K, p_prob,
          obs_ind_v = 1:nrow(x), pred_ind_v = nrow(x) + 1:nrow(x_pred))  
      }
    }
  }
  return(c(m_k_m_draws, list(
    "tilde_p_k_m" = tilde_p_k_m_draws, 
    "mu_tilde_p_k_m" = mu_tilde_p_k_m, 
    "labels" = label_m, 
    "comp_prob" = p_m,
    "post_pred_prob" = post_pred_prob,
    "last_z_k_m" = z_k_m)))
}
