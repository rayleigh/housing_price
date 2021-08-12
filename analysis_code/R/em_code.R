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

general_inverse_logit <- function(tilde_z_k_m) {
  p_k_m <- cbind(exp(tilde_z_k_m), 1)
  p_k_m / rowSums(p_k_m)
}

#Calculates von Mises probability for a list of angles given the mean (m_k)
#and concentration (kappa_k)
#Can be used for a vector y with a single or vector m_k and kappa_k
calc_von_mises_prob <- function(y, m_k, kappa_k, log = F) {
  vm_log_prob = kappa_k * cos(y - m_k) - log(2 * pi * besselI(kappa_k, 0))
  if (log) {
    return(vm_log_prob)
  }
  return(exp(vm_log_prob))
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
          calc_von_mises_prob(y, m_k_m[, k], p_k_m[, k], log = T))) +
      sum(dnorm(tilde_p_k_m[, k], mu_tilde_p_k[k], 0.05, log = T))
  }
  return(log_ll)
}

#Calculate the probability an observation belongs to a cluster for EM on SvM-c
calc_r_k_m <- function(y, p, m_k_m, tilde_z_k_m, tilde_p_k_m, mu_tilde_p_k) {
  p_k_m <- exp(tilde_p_k_m)
  r_k_m <- matrix(0, nrow = nrow(m_k_m), ncol = ncol(m_k_m))
  for (k in 1:ncol(r_k_m)) {
    r_k_m[, k] <- p[k] * 
      calc_von_mises_prob(y, m_k_m[, k], p_k_m[, k])
    # dbeta(y, m_k_m[, k] * p_k_m[, k],
    #       (1 - m_k_m[, k]) * p_k_m[, k])
    # dnorm(tilde_p_k_m[, k], mu_tilde_p_k[k], 0.1) *
    # dnorm(tilde_z_k_m[, k])
  }
  return(r_k_m / rowSums(r_k_m))
}

#Find the z_k_m (GP surface) that maximizes the conditional log likelihood for EM on SvM and SvM-c
#assuming a projected Gaussian Process and a noncentered parametrization
optimize_tilde_z_k_m <- function(tilde_z_k_m, y, r_k_m, p_k_m, mu_z_k, cov_chol) {
  for (k in 1:ncol(r_k_m)) {
    for (j in 2 * (k - 1) + 1:2) {
      next_tilde_z_k <- tilde_z_k_m[, j]
      prev_tilde_z_k <- next_tilde_z_k
      if (j %% 2 == 1) {
        other_z_k = cov_chol %*% tilde_z_k_m[, j + 1] + mu_z_k[j + 1]
      } else {
        other_z_k = cov_chol %*% tilde_z_k_m[, j - 1] + mu_z_k[j - 1]
      }
      i = 0
      repeat {
        next_z_k = cov_chol %*% next_tilde_z_k + mu_z_k[j]
        if (j %% 2 == 1) {
          next_m_k = mapply(angle, next_z_k, other_z_k)
          grad_z_k = -1
        } else {
          next_m_k = mapply(angle, other_z_k, next_z_k)
          grad_z_k = 1
        }
        # grad_z_k <- r_k_m[, k] * next_m_k * (1 - next_m_k) * diag(cov_chol) * (
        #   p_k_m[,k] * log(y) - p_k_m[,k] * log(1 - y) -
        #   p_k_m[,k] * digamma(p_k_m[,k] * next_m_k) + 
        #   p_k_m[,k] * digamma(p_k_m[,k] * (1 - next_m_k))) 
        # grad_z_k <- r_k_m[, k] * next_m_k * (1 - next_m_k) * p_k_m[,k] * (
        #   logit(y) - digamma(p_k_m[,k] * next_m_k) + digamma(p_k_m[,k] * (1 - next_m_k))) 
        grad_z_k <- grad_z_k * r_k_m[, k] * other_z_k / (next_z_k^2 + other_z_k^2) * 
          p_k_m[,k] * sin(y - next_m_k) 
        grad_z_k <- crossprod(cov_chol, grad_z_k)
        # grad_z_k <- grad_z_k - r_k_m[, k] * next_tilde_z_k
        grad_z_k <- grad_z_k - next_tilde_z_k
        next_tilde_z_k = next_tilde_z_k + 1e-4 * grad_z_k
        #print(head(grad_z_k))
        if (norm(next_tilde_z_k - prev_tilde_z_k, type = "2") < 1e-4) {
          break
        }
        if (i > 1000) {
          break
        }
        i = i + 1
        prev_tilde_z_k <- next_tilde_z_k
      }
      tilde_z_k_m[, j] <- next_tilde_z_k 
    }
  }
  return(tilde_z_k_m)
}

#Find the tilde_p_k_m (concentration) that maximizes the conditional log likelihood for EM on SvM and SvM-c
#assuming a hierarchical concentration parameter
optimize_tilde_p_k_m <- function(tilde_p_k_m, y, r_k_m, m_k_m, mu_tilde_p_k) {
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
      grad_tilde_p_k <- r_k_m[, k] * next_p_k * (
        cos(y - m_k_m[, k]) - besselI(next_p_k, -1) / besselI(next_p_k, 0)) 
      # grad_tilde_p_k = grad_tilde_p_k - r_k_m[, k] / (0.1^2) * (next_tilde_p_k - mu_tilde_p_k[k])
      grad_tilde_p_k = grad_tilde_p_k - 1 / (0.05^2) * (next_tilde_p_k - mu_tilde_p_k[k])
      next_tilde_p_k = next_tilde_p_k + 1e-4 * grad_tilde_p_k
      #print(norm(next_tilde_p_k - prev_tilde_p_k, type = "2"))
      if (norm(next_tilde_p_k - prev_tilde_p_k, type = "2") < 1e-4) {
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

#1e-6 for 2 components
#1e-7 for 3 components
#
get_em_hot_start <- function(y, x, mu_z_k, sigma, l, K, mu_tilde_p_k = rep(0, K)) {
  cov <- as.matrix(dist(x))^2
  cov <- sigma^2 * exp(-cov / (2 * l^2)) + diag(rep(1e-9, nrow(cov)))
  cov_chol <- t(chol(cov))
  
  next_p <- rep(1 / K, K)
  prev_p <- next_p
  next_tilde_z_k_m <- matrix(0, nrow = length(y), ncol = K * 2)
  prev_tilde_z_k_m <- next_tilde_z_k_m
  # next_m_k_m <- 2 * pi * invlogit(sweep(cov_chol %*% next_tilde_z_k_m, 2, mu_z_k, "+"))
  next_z_k_m <- sweep(cov_chol %*% next_tilde_z_k_m, 2, mu_z_k, "+")
  next_m_k_m <- sapply(1:K, function(k) {
    mapply(angle, next_z_k_m[, 2 * (k - 1) + 1], next_z_k_m[, 2 * (k - 1) + 2])
  })
  prev_m_k_m <- next_m_k_m
  next_tilde_p_k_m <- matrix(mu_tilde_p_k, nrow = length(y), ncol = K, byrow = T)
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
    r_k_m <- matrix(1, nrow = length(y), ncol = 1)
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
    next_tilde_z_k_m <- optimize_tilde_z_k_m(next_tilde_z_k_m, y, r_k_m, 
                                             exp(next_tilde_p_k_m), mu_z_k, cov_chol)
    next_z_k_m <- sweep(cov_chol %*% next_tilde_z_k_m, 2, mu_z_k, "+")
    next_m_k_m <- sapply(1:K, function(k) {
      mapply(angle, next_z_k_m[, 2 * (k - 1) + 1], next_z_k_m[, 2 * (k - 1) + 2])
    })
    # print("Update tilde_z_k_m")
    
    next_tilde_p_k_m <- optimize_tilde_p_k_m(next_tilde_p_k_m, y, r_k_m, next_m_k_m, next_mu_tilde_p_k)
    next_mu_tilde_p_k <- update_mu_tilde_p_k(next_tilde_p_k_m, r_k_m)
    # print("Update tilde_p_k_m")
    
    #E step
    if (K > 1) {
      r_k_m <- calc_r_k_m(y, next_p, next_m_k_m, next_tilde_z_k_m, next_tilde_p_k_m, next_mu_tilde_p_k)
    } else {
      r_k_m <- matrix(1, nrow = length(y), ncol = 1)
    }
    next_ll = calc_expected_cond_log_likelihood(
      y, r_k_m, next_p, next_m_k_m, next_tilde_z_k_m, cov_chol, 
      next_tilde_p_k_m, next_mu_tilde_p_k)
    # print(prev_ll)
    # print(next_ll)
    if ((next_ll - prev_ll) / abs(prev_ll) < 1e-7) {
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

calc_expected_cond_log_likelihood_gp_mixture_prob <- function(
  y, r_k_m, p_k_m, tilde_z_k_m, m_list, kappa_list) {
  log_ll = sum(dnorm(tilde_z_k_m, log = T)) + 
    sum(dgamma(kappa_list, 1, 1, log = T)) +
    sum(calc_von_mises_prob(m_list, pi, 0, log = T))
  for (k in 1:ncol(p_k_m)) {
    log_ll = log_ll + sum(
      r_k_m[, k] * (
        log(p_k_m[, k]) +
          calc_von_mises_prob(y, m_list[k], kappa_list[k], log = T)))
  }
  return(log_ll)
}

calc_r_k_m_gp_mixture_prob <- function(y, p_k_m, m_list, kappa_list) {
  r_k_m <- matrix(0, nrow = nrow(p_k_m), ncol = ncol(p_k_m))
  for (k in 1:ncol(r_k_m)) {
    r_k_m[, k] <- p_k_m[, k] * 
      calc_von_mises_prob(y, m_list[k], kappa_list[k])
    # dbeta(y, m_k_m[, k] * p_k_m[, k],
    #       (1 - m_k_m[, k]) * p_k_m[, k])
    # dnorm(tilde_p_k_m[, k], mu_tilde_p_k[k], 0.1) *
    # dnorm(tilde_z_k_m[, k])
  }
  return(r_k_m / rowSums(r_k_m))
}

update_tilde_z_k_gp_mixture_prob <- function(prev_tilde_z_k, r_k_m, cov_chol) {
  next_tilde_z_k <- prev_tilde_z_k
  i = 0
  repeat {
    next_p_k_m <- invlogit(cov_chol %*% next_tilde_z_k)
    grad_z_k <- r_k_m[, 1]  - next_p_k_m 
    grad_z_k <- crossprod(cov_chol, grad_z_k)
    # grad_z_k <- grad_z_k - r_k_m[, k] * next_tilde_z_k
    grad_z_k <- grad_z_k - next_tilde_z_k
    next_tilde_z_k = next_tilde_z_k + 1e-3 * grad_z_k
    #print(head(grad_z_k))
    if (norm(next_tilde_z_k - prev_tilde_z_k, type = "2") < 1e-4) {
      break
    }
    if (i > 1000) {
      break
    }
    i = i + 1
    prev_tilde_z_k <- next_tilde_z_k
  }
  return(next_tilde_z_k)
}

update_tilde_z_k_gp_mixture_prob_K <- function(tilde_z_k, r_k_m, cov_chol) {
  for (k in 1:ncol(tilde_z_k)) {
    next_tilde_z_k <- tilde_z_k[, k]
    prev_tilde_z_k <- next_tilde_z_k
    i = 0
    repeat {
      tilde_z_k[, k] <- next_tilde_z_k
      next_p_k_m <- general_inverse_logit(tilde_z_k)
      grad_z_k <- (r_k_m[, k]  - next_p_k_m[, k])
      grad_z_k <- crossprod(cov_chol, grad_z_k)
      # grad_z_k <- grad_z_k - r_k_m[, k] * next_tilde_z_k
      grad_z_k <- grad_z_k - next_tilde_z_k
      next_tilde_z_k = next_tilde_z_k + 1e-4 * grad_z_k
      #print(head(grad_z_k))
      if (norm(next_tilde_z_k - prev_tilde_z_k, type = "2") < 1e-5) {
        break
      }
      if (i > 1000) {
        break
      }
      i = i + 1
      prev_tilde_z_k <- next_tilde_z_k
    }
    tilde_z_k[, k] <- next_tilde_z_k
  }
  return(tilde_z_k)
}


update_kappa_list_gp_mixture_prob <- function(prev_kappa_list, y, r_k_m, m_list, prior = T) {
  next_kappa_list <- prev_kappa_list
  i = 0
  for (k in 1:ncol(r_k_m)) {
    next_kappa = prev_kappa_list[k]
    prev_kappa = next_kappa
    repeat {
      grad_kappa_k <- sum(r_k_m[, k] * 
                            (cos(y - m_list[k]) - besselI(next_kappa, -1) / besselI(next_kappa, 0))) 
      if (prior) {
        grad_kappa_k = grad_kappa_k - 1
      }
      next_kappa = next_kappa + 1e-3 * grad_kappa_k
      #print(head(grad_z_k))
      if (next_kappa < 0) {
        next_kappa = prev_kappa
        break
      }
      if (norm(next_kappa_list[k] - prev_kappa_list[k], type = "2") < 1e-4) {
        break
      }
      if (i > 1000) {
        break
      }
      i = i + 1
      prev_kappa = next_kappa
    }
    next_kappa_list[k] <- next_kappa
  }
  return(next_kappa_list)
}

sort_results <- function(p_k_m, tilde_z_k_m, m_list, kappa_list, ll) {
  sorted_order <- order(m_list)
  tilde_z_k_m <- cbind(tilde_z_k_m, 0)
  tilde_z_k_m <- tilde_z_k_m[, sorted_order]
  tilde_z_k_m <- tilde_z_k_m - tilde_z_k_m[, ncol(tilde_z_k_m)]
  p_k_m <- p_k_m[, sorted_order]
  m_list <- m_list[sorted_order]
  kappa_list <- kappa_list[sorted_order]
  return(list(p_k_m, tilde_z_k_m[, 1:(ncol(p_k_m) - 1)], m_list, kappa_list, ll))
}

get_em_hot_start_gp_mixture_prob <- function(y, x, sigma, l, K) {
  cov <- as.matrix(dist(x))^2
  cov <- sigma^2 * exp(-cov / (2 * l^2)) + diag(rep(1e-9, nrow(cov)))
  cov_chol <- t(chol(cov))
  
  next_tilde_z_k_m <- matrix(0, nrow = length(y), ncol = K - 1)
  prev_tilde_z_k_m <- next_tilde_z_k_m
  next_p_k_m <- general_inverse_logit(cov_chol %*% next_tilde_z_k_m)
  prev_p_k_m <- next_p_k_m
  #p_k_m <- cbind(p_k_m, 1 - p_k_m)
  next_m_list <- 2 * pi / K * (0:(K - 1)) + pi / K 
  prev_m_list <- next_m_list
  next_kappa_list <- rep(1, K)
  prev_kappa_list <- next_kappa_list
  # for (k in 1:K) {
  #   tilde_p_k_m[, k] <- rnorm(length(y), mu_tilde_p_k[k], 0.1)
  # }
  r_k_m <- calc_r_k_m_gp_mixture_prob(y, next_p_k_m, next_m_list, next_kappa_list)
  next_ll = calc_expected_cond_log_likelihood_gp_mixture_prob(
    y, r_k_m, next_p_k_m, next_tilde_z_k_m, next_m_list, next_kappa_list)
  prev_ll = next_ll
  repeat {
    #M step
    if (K == 2) {
      next_tilde_z_k_m <- update_tilde_z_k_gp_mixture_prob(prev_tilde_z_k_m, r_k_m, cov_chol)
    } else {
      next_tilde_z_k_m <- update_tilde_z_k_gp_mixture_prob_K(prev_tilde_z_k_m, r_k_m, cov_chol)
    }
    # p_k_m <- exp(-(cov_chol %*% tilde_z_k_m))
    # p_k_m <- cbind(p_k_m, 1)
    # p_k_m <- p_k_m / rowSums(p_k_m)
    next_p_k_m <- general_inverse_logit(next_tilde_z_k_m)
    # print("Update tilde_z_k_m")
    
    for (k in 1:ncol(r_k_m)) {
      next_m_list[k] = angle(sum(r_k_m[, k] * cos(y)), sum(r_k_m[, k] * sin(y)))
    }
    # if (K > 2) {
    #   sorted_order <- order(next_m_list)
    #   next_tilde_z_k_m <- 
    #   next_m_list <- next_m_list[sorted_order]
    #   next_kappa_list <- next_kappa_list[sorted_order]
    # } else {
    #   next_m_list <- sort(next_m_list)
    # }
    
    
    next_kappa_list <- update_kappa_list_gp_mixture_prob(next_kappa_list, y, r_k_m, next_m_list)
    # print("Update tilde_p_k_m")
    
    #E step
    r_k_m <- calc_r_k_m_gp_mixture_prob(y, next_p_k_m, next_m_list, next_kappa_list)
    next_ll = calc_expected_cond_log_likelihood_gp_mixture_prob(
      y, r_k_m, next_p_k_m, next_tilde_z_k_m, next_m_list, next_kappa_list)
    print(prev_ll)
    print(next_ll)
    if ((next_ll - prev_ll) / abs(prev_ll) < 1e-4) {
      break
    }
    prev_ll = next_ll
    prev_tilde_z_k_m <- next_tilde_z_k_m
    prev_p_k_m <- next_p_k_m
    #p_k_m <- cbind(p_k_m, 1 - p_k_m)
    prev_m_list <- next_m_list
    prev_kappa_list <- next_kappa_list
  }
  if (next_ll > prev_ll) {
    return(sort_results(next_p_k_m, next_tilde_z_k_m, next_m_list, next_kappa_list, next_ll))
  }
  return(sort_results(prev_p_k_m, prev_tilde_z_k_m, prev_m_list, prev_kappa_list, prev_ll))
}

calc_r_k_m_ind_mixture <- function(y, p_k_m, m_list, kappa_list) {
  r_k_m <- sapply(1:length(p_k_m), function(k) {
    p_k_m[k] * calc_von_mises_prob(y, m_list[k], kappa_list[k])
  })
  r_k_m / rowSums(r_k_m)
}

calc_expected_cond_log_likelihood_ind_mixture <- function(
  y, r_k_m, p_k_m, m_list, kappa_list) {
  log_prob = sum(calc_von_mises_prob(m_list, pi, 0, log = T))
  log_prob = log_prob + log(ddirichlet(p_k_m, rep(1, length(p_k_m))))
  for (k in 1:ncol(r_k_m)) {
    log_prob = log_prob + sum(r_k_m[, k] * (log(p_k_m[k]) +
                                              calc_von_mises_prob(y, m_list[k], kappa_list[k], log = T)))
  }
  return(log_prob)
}

get_em_hot_start_ind_mixture <- function(y, x, K) {
  next_p_k_m <- rep(1/K, K)
  prev_p_k_m <- next_p_k_m
  #p_k_m <- cbind(p_k_m, 1 - p_k_m)
  next_m_list <- 2 * pi / K * (0:(K - 1)) + pi / K 
  # next_m_list <- 2 * pi * runif(K)
  # next_m_list <- sort(next_m_list)
  prev_m_list <- next_m_list
  next_kappa_list <- rep(1, K)
  prev_kappa_list <- next_kappa_list
  # for (k in 1:K) {
  #   tilde_p_k_m[, k] <- rnorm(length(y), mu_tilde_p_k[k], 0.1)
  # }
  r_k_m <- calc_r_k_m_ind_mixture(y, next_p_k_m, next_m_list, next_kappa_list)
  next_ll = calc_expected_cond_log_likelihood_ind_mixture(
    y, r_k_m, next_p_k_m, next_m_list, next_kappa_list)
  prev_ll = next_ll
  repeat {
    #M step
    next_p_k_m <- colMeans(r_k_m)
    # print(next_p_k_m)
    # print("Update tilde_z_k_m")
    
    for (k in 1:ncol(r_k_m)) {
      next_m_list[k] = angle(sum(r_k_m[, k] * cos(y)), sum(r_k_m[, k] * sin(y)))
    }
    sorted_order <- order(next_m_list)
    r_k_m <- r_k_m[, sorted_order]
    next_p_k_m <- next_p_k_m[sorted_order]
    next_m_list <- next_m_list[sorted_order]
    next_kappa_list <- next_kappa_list[sorted_order]
    # if (K > 2) {
    #   sorted_order <- order(next_m_list)
    #   next_tilde_z_k_m <- 
    #   next_m_list <- next_m_list[sorted_order]
    #   next_kappa_list <- next_kappa_list[sorted_order]
    # } else {
    #   next_m_list <- sort(next_m_list)
    # }
    next_kappa_list <- update_kappa_list_gp_mixture_prob(
      next_kappa_list, y, r_k_m, next_m_list, prior = F)
    # print("Update tilde_p_k_m")
    
    #E step
    r_k_m <- calc_r_k_m_ind_mixture(y, next_p_k_m, next_m_list, next_kappa_list)
    next_ll = calc_expected_cond_log_likelihood_ind_mixture(
      y, r_k_m, next_p_k_m, next_m_list, next_kappa_list)
    #print(prev_ll)
    #print(next_ll)
    if ((next_ll - prev_ll) / abs(prev_ll) < 1e-4) {
      break
    }
    prev_ll = next_ll
    prev_p_k_m <- next_p_k_m
    #p_k_m <- cbind(p_k_m, 1 - p_k_m)
    prev_m_list <- next_m_list
    prev_kappa_list <- next_kappa_list
  }
  if (next_ll > prev_ll) {
    sorted_order <- order(next_m_list)
    return(list(next_p_k_m[sorted_order], next_m_list[sorted_order], 
                next_kappa_list[sorted_order], next_ll))
  }
  sorted_order <- order(prev_m_list)
  return(list(prev_p_k_m[sorted_order], prev_m_list[sorted_order], 
              prev_kappa_list[sorted_order], prev_ll))
}

