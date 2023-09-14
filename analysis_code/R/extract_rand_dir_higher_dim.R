library(tidyverse)
library(rstan)

transform_from_spherical_coordinates <- function(phi_list) {
  cart_coord <- rep(1, length(phi_list) + 1)
  for (i in 1:length(phi_list)) {
    cart_coord[i:length(cart_coord)] <- cart_coord[i:length(cart_coord)] *
      c(cos(phi_list[i]), rep(sin(phi_list[i]), length(phi_list) - i + 1))
  }
  cart_coord[abs(cart_coord) < .Machine$double.neg.eps] = 0
  return(cart_coord)
}

#Given a point in the positive orthant, converts to spherical cooordinates
get_phi_list <- function(start_row) {
  d = length(start_row)
  norm_sum <- sqrt(rev(cumsum(rev(start_row^2))[-1]))
  phi_list <- acos((start_row[-d] / norm_sum))
  zero_inds <- which(abs(norm_sum) < .Machine$double.eps)
  phi_list[zero_inds] = 0
  if (start_row[d] < -.Machine$double.eps) {
    phi_list[length(phi_list)] = 2 * pi - phi_list[length(phi_list)]
  }
  return(phi_list)
}

build_rotation_m <- function(simplex_row) {
  start_row <- sqrt(simplex_row)
  start_phi <- get_phi_list(start_row)
  return(unname(cbind(
    start_row,
    sapply(1:length(start_phi), function(i) {
      tmp <- start_phi
      if (i > 1) {
        tmp[1:(i - 1)] = pi / 2
      }
      tmp[i] <- tmp[i] + pi / 2
      transform_from_spherical_coordinates(tmp)
    }))))
}

extract_angle_clean_test <- function(
  clean_test, source_year, dest_year, inc_cat) {
 
  #Set up data
  source_info <- clean_test[clean_test$year == source_year, c(inc_cat, "trtid10")]
  dest_info <- clean_test[clean_test$year == dest_year, c(inc_cat, "trtid10")]
  common_tract_ids <- intersect(source_info$trtid10, dest_info$trtid10)
  common_tract_ids <- common_tract_ids[!is.na(common_tract_ids)]
  source_info <- source_info[source_info$trtid10 %in% common_tract_ids,]
  source_info <- source_info %>% arrange(trtid10) %>% data.matrix()
  dest_info <- dest_info[dest_info$trtid10 %in% common_tract_ids,]
  dest_info <- dest_info %>% arrange(trtid10) %>% data.matrix()
  # plot_info <- cbind(source_info[, inc_cat], dest_info[, inc_cat])
  # labels1 <- sapply(inc_cat, function(cat) {paste(cat, 1, sep = "_")})
  # labels2 <- sapply(inc_cat, function(cat) {paste(cat, 2, sep = "_")})
  # colnames(plot_info) <- c(labels1, labels2)
  
  all_observed_norm_dir <- matrix(0, nrow = nrow(source_info),
                                  ncol = length(inc_cat))
  angle_observed_dir <- matrix(0, nrow = nrow(source_info),
                              ncol = length(inc_cat) - 2)  
  observed_dir <- matrix(0, nrow = nrow(source_info),
                         ncol = length(inc_cat) - 1)
  for (i in 1:nrow(source_info)) {
    all_observed_norm_dir[i,] <- 
      solve(build_rotation_m(source_info[i, inc_cat]),
            sqrt(dest_info[i, inc_cat]))
    angle_observed_dir[i, ] <- get_phi_list((all_observed_norm_dir[i,]))[-1]
    observed_dir[i,] <- transform_from_spherical_coordinates(angle_observed_dir[i,])
  }
  return(list(source_info, dest_info, 
              all_observed_norm_dir, angle_observed_dir, 
              observed_dir))
}

gen_geodesic_higher_dim <- lapply(1990:2009, function(year) {
  extract_angle_clean_test(clean_test, year, year + 1, split_names)
})
file_prefix <- "gp_dir_info_von_mises_fisher_no_duplicates_"
for (i in 1:20) {
  m_l = 0.5
  m_s = 0.5
  d = 6
  observed_angle_info <- gen_geodesic_higher_dim[[i]]
  keep_ind = which(abs(1 - observed_angle_info[[3]][,1]) > .Machine$double.eps)
  observed_dir <- observed_angle_info[[5]][keep_ind,]
  simplex_location_m = observed_angle_info[[1]][keep_ind, 1:d]
  census_tract_id = observed_angle_info[[1]][keep_ind, d + 1]
  
  data_dir <- cbind(simplex_location_m, observed_dir)
  duplicated_inds <- which(duplicated(data_dir))
  if (length(duplicated_inds) > 0) {
    data_dir <- data_dir[-duplicated_inds,]
    mod_length = c(mod_length, nrow(data_dir))
  }
  print(nrow(data_dir))
  
  data_pred_ind <- (sample(nrow(data_dir)))[1:(floor(.1 * nrow(data_dir)))]
  simplex_location_m <- data_dir[-data_pred_ind, 1:ncol(simplex_location_m)]
  observed_dir <- data_dir[
    -data_pred_ind, ncol(simplex_location_m) + 1:ncol(observed_dir)]
  N <- length(observed_dir)
  data_census_tract_id <- census_tract_id[-data_pred_ind]
  
  simplex_prediction_m <- data_dir[data_pred_ind, 1:ncol(simplex_location_m)]
  pred_dir <- data_dir[
    data_pred_ind, ncol(simplex_location_m) + 1:ncol(observed_dir)]
  N_pred = nrow(simplex_prediction_m)
  pred_census_tract_id <- census_tract_id[data_pred_ind]
  
  stan_rdump(c("N", "N_pred", "d", "observed_dir",
               "simplex_location_m", "simplex_prediction_m", 
               "pred_dir", "data_census_tract_id", "pred_census_tract_id"),
             file = paste(file_prefix, i, ".R", sep = ""))
}

