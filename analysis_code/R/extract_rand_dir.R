library(tidyverse)
library(rstan)

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

simplex_hellinger_distance <- function(x, y) {
  1 / sqrt(2) * sqrt(sum((sqrt(x) - sqrt(y))^2))
}

build_rotation_matrix_3d <- function(start_row) {
  start_row <- sqrt(start_row)
  theta0 = acos(start_row[3])
  psi0 = atan(start_row[2] / start_row[1])
  if (is.nan(psi0)) {
    return(diag(3))
  }
  sin_psi0 = sin(psi0)
  cos_psi0 = cos(psi0)
  sin_theta0 = sin(theta0)
  cos_theta0 = cos(theta0)
  rotation_m <- cbind(
    c(cos_theta0 * cos_psi0,
      cos_theta0 * sin_psi0,
      -sin_theta0),
    c(-sin_psi0, cos_psi0, 0),
    c(sin_theta0 * cos_psi0,
      sin_theta0 * sin_psi0,
      cos_theta0))
  apply(rotation_m, 2, function(col) {
    col[abs(col) < .Machine$double.eps] = 0
    sign(col) * sqrt(col^2 / sum(col^2))
  })
}

calculate_angle_on_orthant <- function(start_v, end_v) {
  rotation_m <- build_rotation_matrix_3d(start_v)
  reoriented_coords <- solve(rotation_m, sqrt(end_v))
  return(angle(reoriented_coords[1], reoriented_coords[2]))
}

get_hellinger_geodesic_simplex_corner_endpoint_3d <- 
  function(psi1, theta1, psi2, theta2) {
    if (is.nan(psi2)) {
      return(c(0, 0, 1))
    }
    if (is.nan(psi1)) {
      return(c(cos(psi2), sin(psi2), 0))
    }
    if (abs(psi1 - psi2) < .Machine$double.eps) {
      if (theta1 > theta2) {
        return(c(0,0,1))
      }
      return(c(cos(psi1), sin(psi1), 0))
    }
    # if (abs(theta1 - theta2) < .Machine$double.eps) {
    #   if (psi1 > psi2) {
    #     return(c(1,0,0))
    #   }
    #   return(c(0,1,0))
    # }
    # if (psi2 == 0 || psi2 == pi / 2 || theta2 == pi / 2) {
    #   return(c(sin(theta2) * cos(psi2),
    #            sin(theta2) * sin(psi2),
    #            cos(theta2)))
    # }
    return(rep(-1, 3))
  }

#Finds the hellinger geodesic simplex endpoint in the direction of y_row from x_row
#Need to indicate whether the rows are in the simplex or already on the orthant
#Returns an endpoint in the simplex
get_hellinger_geodesic_simplex_endpoint_3d <- function(x_row, y_row, 
                                                       simplex_val = T, rotated = F) {
  if (simplex_val) {
    x_row <- sqrt(x_row)
    y_row <- sqrt(y_row)
  }
  psi0 <- atan((x_row[3] * y_row[1] - y_row[3] * x_row[1]) / (
    y_row[3] * x_row[2] - x_row[3] * y_row[2]))
  if (is.nan(psi0)) {
    if (rotated) {
      return(rep(-1, 3))
    } else {
      trans_x_row <- x_row[c(3, 2, 1)]
      trans_y_row <- y_row[c(3, 2, 1)]
      return((get_hellinger_geodesic_simplex_endpoint_3d(
        trans_x_row, trans_y_row, simplex_val = F, rotated = T))[c(3, 2, 1)])
    }
  }
  if (abs(psi0) < .Machine$double.eps) {
    psi0 = 0
  }
  a = x_row[3] / (x_row[1] * cos(psi0) + x_row[2] * sin(psi0))
  # print(a)
  # print(x_row[1] * cos(psi0) + x_row[2] * sin(psi0))
  # print(y_row[1] * cos(psi0) + y_row[2] * sin(psi0))
  # print(psi0)
  if (is.nan(a) || a == 0) {
    a = y_row[3] / (y_row[1] * cos(psi0) + y_row[2] * sin(psi0))
  }
  theta1 = acos(x_row[3])
  theta2 = acos(y_row[3])
  psi1 = angle(x_row[1], x_row[2])
  psi2 = angle(y_row[1], y_row[2])
  dest_vector <- get_hellinger_geodesic_simplex_corner_endpoint_3d(
    psi1, theta1, psi2, theta2)
  if (all(dest_vector == -1)) {
    final_psi = c(0, pi/2)
    final_theta = atan(1 / (a * cos(final_psi - psi0)))
    if (psi0 < 0 && ((pi/2 + psi0) >= 0)) {
      final_psi = c(final_psi, pi/2 + psi0)
      final_theta = c(final_theta, pi / 2)
    }
    valid_ind <- which(final_theta >= 0 &
                         final_theta <= pi / 2
                       & final_psi != psi1)
    final_psi = final_psi[valid_ind]
    final_theta = final_theta[valid_ind]
    if (length(valid_ind) == 2) {
      if (psi2 > psi1) {
        final_theta = final_theta[which.max(final_psi)]
        final_psi = max(final_psi)
      } else {
        final_theta = final_theta[which.min(final_psi)]
        final_psi = min(final_psi)
      }
    }
    dest_vector <- c(sin(final_theta) * cos(final_psi), 
                     sin(final_theta) * sin(final_psi),
                     cos(final_theta))
  }
  dest_vector = dest_vector^2
  dest_vector[dest_vector <= .Machine$double.eps] = 0
  dest_vector = dest_vector / sum(dest_vector)
  return(dest_vector)
}

get_max_hellinger_geodesic_distance_3d <- function(x_row, y_row) {
  if (all(abs(x_row - y_row) < .Machine$double.eps)) {
    return(0)
  }
  dest_v <- get_hellinger_geodesic_simplex_endpoint_3d(x_row, y_row)
  return(simplex_hellinger_distance(x_row, dest_v))
}


model_check_for_year <- function(clean_test, source_year, dest_year, inc_cat) {
  #Set up data
  source_info <- clean_test[clean_test$year == source_year, c(inc_cat, "trtid10")]
  dest_info <- clean_test[clean_test$year == dest_year, c(inc_cat, "trtid10")]
  common_tract_ids <- intersect(source_info$trtid10, dest_info$trtid10)
  common_tract_ids <- common_tract_ids[!is.na(common_tract_ids)]
  source_info <- source_info[source_info$trtid10 %in% common_tract_ids,]
  source_info <- source_info %>% arrange(trtid10)
  dest_info <- dest_info[dest_info$trtid10 %in% common_tract_ids,]
  dest_info <- dest_info %>% arrange(trtid10)
  plot_info <- cbind(source_info[, inc_cat], dest_info[, inc_cat])
  labels1 <- sapply(inc_cat, function(cat) {paste(cat, 1, sep = "_")})
  labels2 <- sapply(inc_cat, function(cat) {paste(cat, 2, sep = "_")})
  colnames(plot_info) <- c(labels1, labels2)
  
  # #Check uniform rotation
  dist_prop_theta_0 <- rep(0, nrow(source_info))
  dist_prop_beta_0 <- rep(0, nrow(source_info))
  dist_beta_0 <- rep(0, nrow(source_info))
  angle_prop_beta_0 <- rep(0, nrow(source_info))
  max_dist_beta_0 = rep(0, nrow(source_info))
  for (i in 1:nrow(source_info)) {
     beta_0_simplex_dist = simplex_hellinger_distance(
       source_info[i,], dest_info[i,])
     angle_prop_beta_0[i] = calculate_angle_on_orthant(
       source_info[i,], dest_info[i,]) / (2 * pi)
     if (beta_0_simplex_dist > .Machine$double.eps) {
       rotation_m = build_rotation_matrix_3d(source_info[i,])
       dist_prop_theta_0[i] = acos(solve(rotation_m, sqrt(dest_info[i,]))[3]) / 
         acos(solve(rotation_m, sqrt(
           get_hellinger_geodesic_simplex_endpoint_3d(source_info[i,], dest_info[i,])))[3])
       dist_beta_0[i] = beta_0_simplex_dist
       max_dist_beta_0[i] = get_max_hellinger_geodesic_distance_3d(source_info[i,], dest_info[i,])
       dist_prop_beta_0[i] = beta_0_simplex_dist / max_dist_beta_0[i]
     }
  }
  return(list(dist_prop_beta_0, max_dist_beta_0,
              dist_beta_0, angle_prop_beta_0,
              dist_prop_theta_0,
              source_info, dest_info))
              #source_info[which(dist_prop > 1.01),],
              #dest_info[which(dist_prop > 1.01),]))
}

split_names <- c("inc_0_100_c_", "inc_100_200_c", "inc_200_plus_c_")
years <- 1990:2010
geodesic_result <- lapply(2:length(years) - 1, function(i) {
  model_check_for_year(clean_test, years[i], years[i + 1], split_names)
})

file_prefix = "gp_dir_info_von_mises_no_duplicates_"
for (i in 1:20) {
  keep_ind = which(geodesic_result[[i]][[6]] > .Machine$double.eps &
                     (!is.na(geodesic_result[[i]][[6]])))
  observed_dir = geodesic_result[[i]][[4]][keep_ind] * 2 * pi
  observed_dir[observed_dir > (2 * pi - .Machine$double.eps)] = 0
  simplex_location_m = geodesic_result[[i]][[6]][keep_ind,]

  data_dir <- cbind(simplex_location_m, observed_dir)
  duplicated_inds <- which(duplicated(data_dir))
  data_dir <- data_dir[-duplicated_inds,]
  mod_length = c(mod_length, nrow(data_dir))
  print(nrow(data_dir))
  
  data_pred_ind <- (sample(nrow(data_dir)))[1:(floor(.1 * nrow(data_dir)))]
  simplex_location_m <- data_dir[-data_pred_ind, 1:3]
  observed_dir <- data_dir[-data_pred_ind, 4]
  N <- length(observed_dir)
  
  simplex_prediction_m <- data_dir[data_pred_ind, 1:3]
  pred_dir <- data_dir[data_pred_ind, 4]
  N_pred = nrow(simplex_prediction_m)
  d = 3

  stan_rdump(c("N", "N_pred", "d", "observed_dir",
               "simplex_location_m", "simplex_prediction_m", "pred_dir"),
             file = paste("../stan_rdump/gp_dir_info/", file_prefix, i, ".R", sep = ""))
}

