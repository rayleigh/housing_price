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

#Calculate the expected conditional log likelihood for EM on SvM and SvM-c
calc_expected_cond_log_likelihood <- function(
  y, r_k_m, p, m_k_m, tilde_z_k_m, cov_chol, tilde_p_k_m, mu_tilde_p_k) {
  
  p_k_m <- exp(tilde_p_k_m)
  log_ll = sum(dnorm(mu_tilde_p_k, 0, 5, log = T)) +
    sum(dnorm(tilde_z_k_m, log = T))
  for (k in 1:ncol(r_k_m)) {
    # log_ll = log_ll + sum(
    #   r_k_m[, k] * (
    #   log(p[k]) + 
    #     dbeta(y, m_k_m[, k] * p_k_m[, k], 
    #                      (1 - m_k_m[, k]) * p_k_m[, k], log = T) +
    #     dnorm(tilde_p_k_m[, k], mu_tilde_p_k[k], 0.1, log = T) +
    #     dnorm(tilde_z_k_m[, k], log = T)))
    # log_ll = log_ll + sum(
    #   r_k_m[, k] * (
    #     log(p[k]) + 
    #       dbeta(y, m_k_m[, k] * p_k_m[, k], 
    #             (1 - m_k_m[, k]) * p_k_m[, k], log = T))) +
    #       sum(dnorm(tilde_p_k_m[, k], mu_tilde_p_k[k], 0.1, log = T)) 
    log_ll = log_ll + sum(
      r_k_m[, k] * (
        log(p[k]) +
          calc_von_mises_fisher_prob(y, m_k_m[, , k], p_k_m[, k], log = T))) +
      sum(dnorm(tilde_p_k_m[, k], mu_tilde_p_k[k], 0.05, log = T))
  }
  return(log_ll)
}

#Calculate the probability an observation belongs to a cluster for EM on SvM-c
calc_r_k_m <- function(y, p, m_k_m, tilde_z_k_m, tilde_p_k_m, mu_tilde_p_k) {
  p_k_m <- exp(tilde_p_k_m)
  r_k_m <- matrix(0, nrow = nrow(m_k_m), ncol = length(p))
  for (k in 1:ncol(r_k_m)) {
    r_k_m[, k] <- log(p[k]) + 
      calc_von_mises_fisher_prob(y, m_k_m[, , k], p_k_m[, k], log = T)
    # dbeta(y, m_k_m[, k] * p_k_m[, k],
    #       (1 - m_k_m[, k]) * p_k_m[, k])
    # dnorm(tilde_p_k_m[, k], mu_tilde_p_k[k], 0.1) *
    # dnorm(tilde_z_k_m[, k])
  }
  norm_factor <- apply(r_k_m, 1, logSumExp)
  return(exp(r_k_m - norm_factor))
}

#Find the z_k_m (GP surface) that maximizes the conditional log likelihood for EM on SvM and SvM-c
#assuming a projected Gaussian Process and a noncentered parametrization
optimize_tilde_z_k_m_higher_dim <- function(
  tilde_z_k_m, y, r_k_m, p_k_m, mu_z_k, cov_chol, step_size = 1e-4) {
  
  for (k in 1:ncol(r_k_m)) {
    next_tilde_z_k <- tilde_z_k_m[,,k]
    prev_tilde_z_k <- next_tilde_z_k
    i = 0
    repeat {
      next_z_k <- sweep(cov_chol %*% next_tilde_z_k, 2, mu_z_k[k,], "+")
      next_m_k <- next_z_k / apply(next_z_k, 1, function(row) {norm(row, type = "2")})
      grad_z_k <- matrix(0, nrow = nrow(next_m_k), ncol = ncol(next_m_k))
      for (j in 1:nrow(next_z_k)) {
        z_k_norm <- norm(next_z_k, type = "2")
        grad_z_k[j,] <- r_k_m[j, k] * p_k_m[j, k] * 
          y[j,] %*% ((diag(ncol(y))) / z_k_norm - 
                       tcrossprod(next_z_k[j,]) / z_k_norm^3)
      }
      # grad_z_k <- r_k_m[, k] * (p_k_m[,k] * mu) other_z_k / (next_z_k^2 + other_z_k^2) * 
      #   p_k_m[,k] * sin(y - next_m_k) 
      grad_z_k <- crossprod(cov_chol, grad_z_k)
      # grad_z_k <- grad_z_k - r_k_m[, k] * next_tilde_z_k
      grad_z_k <- grad_z_k - next_tilde_z_k
      next_tilde_z_k = next_tilde_z_k + step_size * grad_z_k
      #print(head(grad_z_k))
      if (norm(next_tilde_z_k - prev_tilde_z_k, type = "2") < step_size) {
        break
      }
      if (i > 1000) {
        break
      }
      i = i + 1
      prev_tilde_z_k <- next_tilde_z_k
    }
    tilde_z_k_m[,, k] <- next_tilde_z_k 
  }
  return(tilde_z_k_m)
}

optimize_tilde_p_k_m_higher_dim <- function(
    tilde_p_k_m, y_m, r_k_m, m_k_m_array, mu_tilde_p_k, step_size = 1e-4) {
  
  d = ncol(y_m)
  for (k in 1:ncol(tilde_p_k_m)) {
    next_tilde_p_k <- tilde_p_k_m[,k]
    prev_tilde_p_k <- next_tilde_p_k
    i = 0
    repeat {
      next_p_k = exp(next_tilde_p_k)
      # grad_tilde_p_k <- r_k_m[, k] * next_p_k * (
      #   m_k_m[,k] * log(y) + (1 - m_k_m[,k]) * log(1 - y) -
      #   m_k_m[,k] * digamma(next_p_k * m_k_m[,k]) - 
      #   (1 - m_k_m[,k]) * digamma(next_p_k * (1 - m_k_m[,k])) +
      #   digamma(next_p_k)) 
      
      # ll_grad <- rowSums(m_k_m_array[,,k] * y_m)
      # ll_grad <- ll_grad + (d / 2 - 1) / next_p_k
      # ll_grad <- ll_grad - (d / 2 - 1) / (next_p_k) + 
      #   besselI(next_p_k, d / 2) / besselI(next_p_k, d / 2 - 1)
      ll_grad <- rowSums(m_k_m_array[,,k] * y_m) - 
        besselI(next_p_k, d / 2) / besselI(next_p_k, d / 2 - 1)
      
      # ll_grad <- ll_grad - (
      #   besselI(next_p_k, d / 2 - 1) * (d / 2 - 1) / (next_p_k) + 
      #     besselI(next_p_k, d / 2)) / 
      #   besselI(next_p_k, d / 2 - 1)

      grad_tilde_p_k <- r_k_m[, k] * next_p_k * ll_grad
      # grad_tilde_p_k = grad_tilde_p_k - r_k_m[, k] / (0.1^2) * (next_tilde_p_k - mu_tilde_p_k[k])
      grad_tilde_p_k = grad_tilde_p_k - 1 / (0.05^2) * (next_tilde_p_k - mu_tilde_p_k[k])
      next_tilde_p_k = next_tilde_p_k + step_size * grad_tilde_p_k
      #print(norm(next_tilde_p_k - prev_tilde_p_k, type = "2"))
      if (norm(next_tilde_p_k - prev_tilde_p_k, type = "2") < step_size) {
        break
      }
      if (i > 1000) {
        break
      }
      i = i + 1
      prev_tilde_p_k <- next_tilde_p_k
    }
    tilde_p_k_m[, k] <- next_tilde_p_k
  }
  return(tilde_p_k_m)
}

#Find the mu_tilde_p_k_m (concentration) that maximizes the conditional log likelihood for EM on SvM and SvM-c
#assuming a hierarchical concentration parameter set-up
update_mu_tilde_p_k <- function(tilde_p_k_m, r_k_m) {
  mu_tilde_p_k <- rep(0, ncol(tilde_p_k_m))
  for (k in 1:ncol(tilde_p_k_m)) {
    # mu_tilde_p_k[k] <- sum(tilde_p_k_m[, k] * r_k_m[, k] / 0.1^2) /
    #   (sum(r_k_m[, k]) / 0.1^2 + 1 / 5^2)
    mu_tilde_p_k[k] <- sum(tilde_p_k_m[, k] / 0.05^2) /
      (length(r_k_m[, k]) / 0.05^2 + 1 / 5^2)
  }
  return(mu_tilde_p_k)
}

matern_kernel_three_halves <- function(x, l) {
  cov <- as.matrix(dist(x))
  cov <- sqrt(3) * cov / l
  (1 + cov) * exp(-cov)
}

matern_kernel_five_halves <- function(x, l) {
  cov <- as.matrix(dist(x))
  cov <- sqrt(5) * cov / l
  (1 + cov + 1/3 * cov^2) * exp(-cov)
}

matern_kernel_one_half <- function(x, l) {
  cov <- as.matrix(dist(x))
  exp(-cov / l)
}

#1e-3 for SvM
#1e-4 for SvM-c
get_em_hot_start_higher_dim <- function(y, x, mu_z_k, sigma, l, K, mu_tilde_p_k = rep(0, K),
                             matern_cov = F, v = 5/2, step_size = 1e-4) {
  if (matern_cov) {
    if (v == 5/2)
      cov <- matern_kernel_five_halves(x, l)
    if (v == 3/2)
      cov <- matern_kernel_three_halves(x, l)
  } else {
    cov <- as.matrix(dist(x))^2
    cov <- sigma^2 * exp(-cov / (2 * l^2))
  }
  cov <- cov + diag(rep(1e-9, nrow(cov)))
  cov_chol <- t(chol(cov))
  next_p <- rep(1 / K, K)
  prev_p <- next_p
  next_tilde_z_k_m <- array(0, dim = c(nrow(y), ncol(y), K))
  prev_tilde_z_k_m <- next_tilde_z_k_m
  # next_m_k_m <- 2 * pi * invlogit(sweep(cov_chol %*% next_tilde_z_k_m, 2, mu_z_k, "+"))
  next_z_k_m <- next_tilde_z_k_m
  next_m_k_m <- next_tilde_z_k_m
  for (k in 1:K) {
    next_z_k_m[,,k] <- sweep(cov_chol %*% next_tilde_z_k_m[,,k], 2, mu_z_k[k,], "+")
    next_m_k_m[,,k] <- next_z_k_m[,,k] / 
      apply(next_z_k_m[,,k], 1, norm, type = "2")
  }
  prev_m_k_m <- next_m_k_m
  next_tilde_p_k_m <- matrix(mu_tilde_p_k, nrow = nrow(y), ncol = K, byrow = T)
  prev_tilde_p_k_m <- next_tilde_p_k_m
  
  next_mu_tilde_p_k <- mu_tilde_p_k
  prev_mu_tilde_p_k <- next_mu_tilde_p_k
  
  # for (k in 1:K) {
  #   tilde_p_k_m[, k] <- rnorm(length(y), mu_tilde_p_k[k], 0.1)
  # }
  if (K > 1) {
    r_k_m <- calc_r_k_m(y, next_p, next_m_k_m, next_tilde_z_k_m, 
                        next_tilde_p_k_m, next_mu_tilde_p_k)
  } else {
    r_k_m <- matrix(1, nrow = nrow(y), ncol = 1)
  }
  next_ll = calc_expected_cond_log_likelihood(
    y, r_k_m, next_p, next_m_k_m, next_tilde_z_k_m, cov_chol, 
    next_tilde_p_k_m, next_mu_tilde_p_k)
  prev_ll = next_ll
  repeat {
    #M step
    #print(next_p)
    if (K > 1) {
      next_p <- colMeans(r_k_m)
      #next_p[K] <- 1 - sum(next_p[1:(K - 1)])
      if (any(next_p < 1e-12)) {
        zero_ind <- which(next_p < 1e-12)
        next_p[zero_ind] = 1e-12
        next_p = next_p / sum(next_p)
        break
      }
    }
    
    next_tilde_z_k_m <- optimize_tilde_z_k_m_higher_dim(
      next_tilde_z_k_m, y, r_k_m,  exp(next_tilde_p_k_m), mu_z_k, cov_chol, step_size)
    
    next_z_k_m <- next_tilde_z_k_m
    next_m_k_m <- next_tilde_z_k_m
    for (k in 1:K) {
      next_z_k_m[,,k] <- sweep(cov_chol %*% next_tilde_z_k_m[,,k], 2, mu_z_k[k,], "+")
      next_m_k_m[,,k] <- next_z_k_m[,,k] / 
        apply(next_z_k_m[,,k], 1, norm, type = "2")
    }
    # next_z_k_m <- sweep(cov_chol %*% next_tilde_z_k_m, 2, mu_z_k, "+")
    # next_m_k_m <- sapply(1:K, function(k) {
    #   mapply(angle, next_z_k_m[, 2 * (k - 1) + 1], next_z_k_m[, 2 * (k - 1) + 2])
    # })
    # print("Update tilde_z_k_m")
    
    next_tilde_p_k_m <- optimize_tilde_p_k_m_higher_dim(
      next_tilde_p_k_m, y, r_k_m, next_m_k_m, next_mu_tilde_p_k, step_size)
    
    next_mu_tilde_p_k <- update_mu_tilde_p_k(next_tilde_p_k_m, r_k_m)
    # print("Update tilde_p_k_m")
    
    #E step
    if (K > 1) {
      r_k_m <- calc_r_k_m(
        y, next_p, next_m_k_m, next_tilde_z_k_m, next_tilde_p_k_m, next_mu_tilde_p_k)
    }

    next_ll = calc_expected_cond_log_likelihood(
      y, r_k_m, next_p, next_m_k_m, next_tilde_z_k_m, cov_chol, 
      next_tilde_p_k_m, next_mu_tilde_p_k)
    print("Update")
    print(prev_ll)
    print(next_ll)
    if ((next_ll - prev_ll) / abs(prev_ll) < (step_size * 1e-2)) {
      break
    }
    prev_ll = next_ll
    prev_p <- next_p
    prev_tilde_z_k_m <- next_tilde_z_k_m
    prev_m_k_m <- next_m_k_m
    prev_tilde_p_k_m <- next_tilde_p_k_m
    prev_mu_tilde_p_k <- next_mu_tilde_p_k
  }
  if (next_ll > prev_ll) {
    return(list(next_p, next_m_k_m, next_tilde_z_k_m, next_tilde_p_k_m, 
                next_mu_tilde_p_k, next_ll))
  }
  return(list(prev_p, prev_m_k_m, prev_tilde_z_k_m, prev_tilde_p_k_m, 
              prev_mu_tilde_p_k, prev_ll))
}