# ============================================================================
# Estimating equation functions for the Joint Frailty Model (JFM)
# ============================================================================


#' Counting Process of Recurrent Event
#'
#' Computes the value of the counting process N_i(t), i.e., the number of
#' recurrent events that occurred at or before time t.
#'
#' @param t A scalar time point.
#' @param recur_time A numeric vector of recurrent event times for a subject.
#'
#' @return Integer count of recurrent events at or before time t.
#'
#' @keywords internal
jfm_counting_process <- function(t, recur_time) {
  if (length(recur_time) == 0) {
    return(0)
  }

# rank t and recurrent time
  a <- c(t, recur_time)
  a <- a[order(a)]
  # find the number of events before t
  return(max(which(a == t)) - 1)
}


#' W_i(t) Evaluated at Each Death Time for Each Subject
#'
#' Constructs the weight matrix W_i(t) evaluated at each death event time
#' for each subject, using the pseudo data set approach.
#'
#' @param theta Frailty variance parameter.
#' @param alpha Coefficient vector for the death sub-model.
#' @param t.start Vector of interval start times.
#' @param I Vector of subject indicators for each pseudo entry.
#' @param Z List of covariate matrices, one per subject.
#' @param td Vector of death event times.
#' @param lambda0_d Baseline hazard point masses for death.
#' @param td.id Subject IDs corresponding to each death event.
#'
#' @return A list with components:
#'   \item{wt_matrix}{Matrix of W_i(t) values (rows = subjects, cols = death times).}
#'   \item{td_id}{Reordered subject IDs for death events.}
#'   \item{index_death_matrix}{Matrix of pseudo entry indices (rows = subjects, cols = death times).}
#'   \item{pseudo_entries}{Sorted pseudo data set.}
#'
#' @keywords internal
jfm_wt_death <- function(theta, alpha, t.start, I, Z, td, lambda0_d, td.id) {
  # construct the pseudo data set (Lj,Ij,zj)
  Z1 <- numeric(0)
  for (i in 1:length(Z)) {
    Z1 <- rbind(Z1, Z[[i]])
  }
  pseudo <- cbind(t.start, I, Z1)
  pseudo <- pseudo[order(pseudo[, 1]), ]
  # rank the death events and also reorder the corresponding point mass (retain the rank-index mapping)
  map_td <- order(td) # ranking to index of td
  td_ranked <- td[map_td]
  td.id <- td.id[map_td] # index of td to id

  lambda0_d_ranked <- lambda0_d # lambda0_d_ranked=lambda0_d[map_td]
  # define a vector to store the index of the pseudo entries whose L is the largest L <= the death time
  lar_L <- numeric(length(td))
  j <- 1
  L <- pseudo[1, 1]
  for (i in 1:length(td)) {
    while (L <= td_ranked[i] && j < nrow(pseudo)) {
      j <- j + 1
      L <- pseudo[j, 1]
    }
    if (L <= td_ranked[i] && j == nrow(pseudo)) {
      lar_L[i] <- j
    } else {
      j <- j - 1
      lar_L[i] <- j
    }
    L <- pseudo[j, 1]
  }
  # define a matrix to store the wt for different subject i and death time (column for death time and row for subject)
  wt_matrix <- matrix(NA, nrow = length(Z), ncol = length(td))
  # also define a matrix to store the index of the pseudo entry whose L is the largest L before each death time for each subject (column for death time and row for subject)
  index_death_matrix <- matrix(NA, nrow = length(Z), ncol = length(td))
  for (i in 1:length(Z)) {
    integral <- 0
    for (j in 1:length(td)) {
      k <- lar_L[j]
      while (pseudo[k, 2] != i) {
        k <- k - 1
      }
      index_death_matrix[i, j] <- k
      integral <- integral + lambda0_d_ranked[j] * exp(alpha %*% pseudo[k, 3:ncol(pseudo)])
      wij <- 1 / (1 + theta * integral)
      wt_matrix[i, j] <- wij
    }
  }
  return(list(wt_matrix = wt_matrix, td_id = td.id, index_death_matrix = index_death_matrix, pseudo_entries = pseudo))
}


#' W_i(t) Evaluated at Each Recurrent Event Time for Each Subject
#'
#' Constructs the weight matrix W_i(t) evaluated at each recurrent event time
#' for each subject, by mapping recurrent event times to the nearest death
#' time in the weight matrix.
#'
#' @param tr Vector of recurrent event times.
#' @param wt_matrix Weight matrix from \code{jfm_wt_death}.
#' @param td Vector of death event times.
#'
#' @return Matrix of W_i(t) values at recurrent event times
#'   (rows = subjects, cols = recurrent event times).
#'
#' @keywords internal
jfm_wt_recurrent <- function(tr, wt_matrix, td) {
  td_ranked <- td[order(td)]
  tr_ranked <- tr[order(tr)]
  j <- 1
  # find the index of the largest death event before each recurrent event
  de_index <- numeric(length(tr))
  flag <- numeric(length(tr))
  for (i in 1:length(tr)) {
    flag[i] <- 1
    t <- tr_ranked[i]
    T <- td_ranked[j]
    while (t > T && j < length(td_ranked)) {
      j <- j + 1
      T <- td_ranked[j]
    }
    if (t > T && j == length(td_ranked)) {
      de_index[i] <- j
    } else {
      if (t == T) {
        de_index[i] <- j
        j <- j + 1
      } else {
        if (j == 1) {
          de_index[i] <- j
          flag[i] <- 0
        } else {
          de_index[i] <- j - 1
        }
      }
    }
  }
  # define the wt_recurrent_subject matrix (row represents subject and column represents the recurrent event time)
  wt_recurrent_subject <- matrix(NA, nrow = nrow(wt_matrix), ncol = length(tr))
  for (i in 1:nrow(wt_matrix)) {
    for (j in 1:length(tr)) {
      if (flag[j] == 0) {
        wt_recurrent_subject[i, j] <- 1
      } else {
        wt_recurrent_subject[i, j] <- wt_matrix[i, de_index[j]]
      }
    }
  }
  return(wt_recurrent_subject)
}


#' Integral in R2i(t) Evaluated at Each Recurrent Event Time
#'
#' Computes the integral component of R2i(t) for each subject at each
#' recurrent event time, using the pseudo data set approach.
#'
#' @param t.start Vector of interval start times.
#' @param I Vector of subject indicators for each pseudo entry.
#' @param Z List of covariate matrices, one per subject.
#' @param beta Coefficient vector for the recurrent event sub-model.
#' @param tr Vector of recurrent event times.
#' @param lambda0_r Baseline hazard point masses for recurrence.
#' @param tr.id Subject IDs corresponding to each recurrent event.
#'
#' @return A list with components:
#'   \item{integral_matrix}{Matrix of integral values (rows = subjects, cols = recurrent times).}
#'   \item{index_recurrent_matrix}{Matrix of pseudo entry indices.}
#'   \item{tr_id}{Reordered subject IDs for recurrent events.}
#'
#' @keywords internal
jfm_r2i_integral <- function(t.start, I, Z, beta, tr, lambda0_r, tr.id) {
  # construct the pseudo data set (Lj,Ij,zj)
  Z1 <- numeric(0)
  for (i in 1:length(Z)) {
    Z1 <- rbind(Z1, Z[[i]])
  }
  pseudo <- cbind(t.start, I, Z1)
  pseudo <- pseudo[order(pseudo[, 1]), ]
  # rank the recurrent events and also reorder the corresponding point mass (retain the rank-index mapping)
  map_tr <- order(tr)
  tr_ranked <- tr[map_tr]
  tr_id <- tr.id[map_tr] # rank to id mapping
  lambda0_r_ranked <- lambda0_r # lambda0_r_ranked=lambda0_r[map_tr]
  # define a matrix to store the index of the pseudo entries whose L is the largest L <= the recurrent time for subject i
  lar_L <- matrix(NA, nrow = length(Z), ncol = length(tr))
  j <- 1
  L <- pseudo[1, 1]

  # Create empty vectors to store indices of pseudo entries which belongs to subject j and are before tr_ranked[i]
  vectors <- vector("list", length(Z))
  for (i in 1:length(tr)) {
    while (L < tr_ranked[i]) {
      vectors[[pseudo[j, 2]]] <- c(vectors[[pseudo[j, 2]]], j)
      j <- j + 1
      L <- pseudo[j, 1]
    }


    for (k in 1:length(Z)) {
      lar_L[k, i] <- vectors[[k]][length(vectors[[k]])]
    }
  }
  # define a matrix to store the integral in r2i for different subject i and recurrent time (column for recurrent time and row for subject)
  integral_matrix <- matrix(NA, nrow = length(Z), ncol = length(tr))
  for (i in 1:length(Z)) {
    integral <- 0
    for (j in 1:length(tr)) {
      k <- lar_L[i, j]
      integral <- integral + lambda0_r_ranked[j] * exp(beta %*% pseudo[k, 3:ncol(pseudo)])
      integral_matrix[i, j] <- integral
    }
  }
  return(list(integral_matrix = integral_matrix, index_recurrent_matrix = lar_L, tr_id = tr_id))
}


#' R2i(t) Evaluated at Each Death Time for Each Subject
#'
#' Computes R2i(t) for each subject at each death event time, by combining
#' the recurrent event integral with the weight matrix.
#'
#' @param t.start Vector of interval start times.
#' @param I Vector of subject indicators for each pseudo entry.
#' @param Z List of covariate matrices, one per subject.
#' @param beta Coefficient vector for the recurrent event sub-model.
#' @param tr Vector of recurrent event times.
#' @param lambda0_r Baseline hazard point masses for recurrence.
#' @param td Vector of death event times.
#' @param wt_matrix Weight matrix from \code{jfm_wt_death}.
#' @param tr.id Subject IDs corresponding to each recurrent event.
#'
#' @return A list with components:
#'   \item{r2i_death_subject_matrix}{Matrix of R2i values at death times.}
#'   \item{index_recurrent_matrix}{Matrix of pseudo entry indices.}
#'   \item{tr_id}{Reordered subject IDs for recurrent events.}
#'
#' @keywords internal
jfm_r2i_death <- function(t.start, I, Z, beta, tr, lambda0_r, td, wt_matrix, tr.id) {
  B <- jfm_r2i_integral(t.start, I, Z, beta, tr, lambda0_r, tr.id)
  integral_matrix <- B$integral_matrix
  index_recurrent_matrix <- B$index_recurrent_matrix
  tr_id <- B$tr_id
  td_ranked <- td[order(td)]
  tr_ranked <- tr[order(tr)]
  r2i_death_subject_matrix <- matrix(NA, nrow = nrow(integral_matrix), ncol = length(td))
  recu_index <- numeric(length(td))
  flag <- numeric(length(td))
  j <- 1
  for (k in 1:length(td)) {
    flag[k] <- 1
    death <- td_ranked[k]
    T <- tr_ranked[j]
    while (death > T && j < length(tr_ranked)) {
      j <- j + 1
      T <- tr_ranked[j]
    }
    if (death > T && j == length(tr_ranked)) {
      recu_index[k] <- j
    } else {
      if (death == T) {
        recu_index[k] <- j
        j <- j + 1
      } else {
        if (j == 1) {
          recu_index[k] <- j
          flag[k] <- 0
        } else {
          recu_index[k] <- j - 1
        }
      }
    }
  }

  for (i in 1:nrow(integral_matrix)) {
    for (j in 1:length(td)) {
      r2i_death_subject_matrix[i, j] <- flag[j] * wt_matrix[i, j] * integral_matrix[i, recu_index[j]]
    }
  }
  return(list(r2i_death_subject_matrix = r2i_death_subject_matrix, index_recurrent_matrix = index_recurrent_matrix, tr_id = tr_id))
}


#' At-Risk Indicator Y*(t) with Status-Dependent Comparison
#'
#' Returns 1 if the subject is at risk at time t, using strict or
#' non-strict inequality depending on the status indicator.
#'
#' @param t A scalar time point.
#' @param y The subject's observed time.
#' @param stat The subject's status indicator (1 = event, 0 = censored).
#'
#' @return Integer 0 or 1 indicator.
#'
#' @keywords internal
jfm_yt_star <- function(t, y, stat) {
  if (stat == 1) {
    if (y > t) {
      result <- 1
    } else {
      result <- 0
    }
  } else {
    result <- ifelse((y >= t), 1, 0)
  }

  return(result)
}


#' G-bar(t) Evaluated at Each Death Time
#'
#' Computes the ratio G-bar(t) at each death event time, used in the
#' estimating equation for theta.
#'
#' @param r2i_death_subject_matrix Matrix of R2i values at death times.
#' @param td Vector of death event times.
#' @param Y Vector of composite censoring/death times per subject.
#' @param STATUS Vector of death indicators per subject.
#' @param list_recur List of recurrent event times per subject.
#'
#' @return A matrix (1 row) of G-bar values, one per death time.
#'
#' @keywords internal
jfm_gt_bar_death <- function(r2i_death_subject_matrix, td, Y, STATUS, list_recur) {
  td_ranked <- td[order(td)]
  # define a vector to store the G_bar value for different death event time
  Gt_bar_death_matrix <- matrix(NA, nrow = 1, ncol = ncol(r2i_death_subject_matrix))
  n <- nrow(r2i_death_subject_matrix)
  for (k in 1:length(td_ranked)) {
    if (r2i_death_subject_matrix[1, k] == 0) {
      Gt_bar_death_matrix[k] <- 0
    } else {
      num <- 0
      denom <- 0
      for (i in 1:n) {
        num <- num + (1 / r2i_death_subject_matrix[i, k]) * jfm_yt_star(td_ranked[k], Y[i], STATUS[i]) * jfm_counting_process(td_ranked[k], list_recur[[i]])
        denom <- denom + jfm_yt_star(td_ranked[k], Y[i], STATUS[i])
      }
      Gt_bar_death_matrix[k] <- num / denom
    }
  }
  return(Gt_bar_death_matrix)
}


#' Estimating Equation for Theta
#'
#' Evaluates the estimating equation for the frailty variance parameter theta.
#'
#' @param td_id Reordered subject IDs for death events.
#' @param r2i_death_subject_matrix Matrix of R2i values at death times.
#' @param td Vector of death event times.
#' @param Y Vector of composite censoring/death times per subject.
#' @param STATUS Vector of death indicators per subject.
#' @param list_recur List of recurrent event times per subject.
#' @param theta Current frailty variance parameter.
#' @param num_recur Integer vector of recurrent event counts per subject.
#'
#' @return Scalar value of the estimating equation for theta.
#'
#' @keywords internal
jfm_est_theta <- function(td_id, r2i_death_subject_matrix, td, Y, STATUS, list_recur, theta, num_recur) {
  Gt_bar_death_matrix <- jfm_gt_bar_death(r2i_death_subject_matrix, td, Y, STATUS, list_recur)
  est <- 0
  for (i in 1:length(td)) {
    est <- est + num_recur[td_id[i]] - (theta + 1) * Gt_bar_death_matrix[i] * r2i_death_subject_matrix[td_id[i], i]
  }
  return(est)
}


#' At-Risk Indicator Y(t)
#'
#' Returns 1 if the subject's observed time y is at least t, 0 otherwise.
#'
#' @param y The subject's observed time.
#' @param t A scalar time point.
#'
#' @return Integer 0 or 1 indicator.
#'
#' @keywords internal
jfm_yt <- function(y, t) {
  if (y >= t) {
    a <- 1
  } else {
    a <- 0
  }
  return(a)
}


#' S1(t) for Death Sub-Model
#'
#' Computes the weighted first moment S1(t) at each death event time for
#' the death sub-model.
#'
#' @param Y Vector of composite censoring/death times per subject.
#' @param wt_matrix Weight matrix from \code{jfm_wt_death}.
#' @param td Vector of death event times.
#' @param index_death_matrix Matrix of pseudo entry indices for death times.
#' @param pseudo_entries Sorted pseudo data set.
#' @param alpha Coefficient vector for the death sub-model.
#'
#' @return Matrix of S1(t) values (rows = covariates, cols = death times).
#'
#' @keywords internal
jfm_s1t_death <- function(Y, wt_matrix, td, index_death_matrix, pseudo_entries, alpha) {
  td_ranked <- td[order(td)]
  nc <- ncol(pseudo_entries)
  # define a vector to store S1t value for different death time
  S1t_death <- matrix(NA, nrow = nc - 2, ncol = length(td))
  for (i in 1:length(td)) {
    sum <- 0
    for (j in 1:length(Y)) {
      sum <- sum + jfm_yt(Y[j], td_ranked[i]) * wt_matrix[j, i] * pseudo_entries[index_death_matrix[j, i], 3:nc] * drop(exp(alpha %*% pseudo_entries[index_death_matrix[j, i], 3:nc]))
    }
    S1t_death[, i] <- sum
  }
  return(S1t_death / length(Y))
}


#' S1(t) for Death Sub-Model with Cross-Fitting
#'
#' Computes the weighted first moment S1(t) at each death event time for
#' the death sub-model, using fold-specific alpha coefficients for cross-fitting.
#'
#' @param Y Vector of composite censoring/death times per subject.
#' @param wt_matrix Weight matrix from \code{jfm_wt_death}.
#' @param td Vector of death event times.
#' @param index_death_matrix Matrix of pseudo entry indices for death times.
#' @param pseudo_entries Sorted pseudo data set.
#' @param alpha_mat Matrix of alpha coefficients (rows = folds).
#' @param CV_map Two-column matrix mapping fold index to subject index.
#'
#' @return Matrix of S1(t) values (rows = covariates, cols = death times).
#'
#' @keywords internal
jfm_s1t_death_cf <- function(Y, wt_matrix, td, index_death_matrix, pseudo_entries, alpha_mat, CV_map) {
  td_ranked <- td[order(td)]
  nc <- ncol(pseudo_entries)
  # define a vector to store S1t value for different death time
  S1t_death <- matrix(NA, nrow = nc - 2, ncol = length(td))
  for (i in 1:length(td)) {
    sum <- 0
    for (j in 1:length(Y)) {
      fold_j <- CV_map[CV_map[, 2] == j, 1]
      alpha <- alpha_mat[fold_j, ]
      sum <- sum + jfm_yt(Y[j], td_ranked[i]) * wt_matrix[j, i] * pseudo_entries[index_death_matrix[j, i], 3:nc] * drop(exp(alpha %*% pseudo_entries[index_death_matrix[j, i], 3:nc]))
    }
    S1t_death[, i] <- sum
  }
  return(S1t_death / length(Y))
}


#' S1(t) for Recurrent Event Sub-Model
#'
#' Computes the weighted first moment S1(t) at each recurrent event time
#' for the recurrence sub-model.
#'
#' @param Y Vector of composite censoring/death times per subject.
#' @param wt_recurrent_subject Weight matrix at recurrent event times.
#' @param tr Vector of recurrent event times.
#' @param index_recurrent_matrix Matrix of pseudo entry indices for recurrent times.
#' @param pseudo_entries Sorted pseudo data set.
#' @param beta Coefficient vector for the recurrent event sub-model.
#'
#' @return Matrix of S1(t) values (rows = covariates, cols = recurrent times).
#'
#' @keywords internal
jfm_s1t_recurrent <- function(Y, wt_recurrent_subject, tr, index_recurrent_matrix, pseudo_entries, beta) {
  tr_ranked <- tr[order(tr)]
  nc <- ncol(pseudo_entries)
  # define a vector to store S1t value for different death time
  S1t_recurrent <- matrix(NA, nrow = nc - 2, ncol = length(tr))
  for (i in 1:length(tr)) {
    sum <- 0
    for (j in 1:length(Y)) {
      sum <- sum + jfm_yt(Y[j], tr_ranked[i]) * wt_recurrent_subject[j, i] * pseudo_entries[index_recurrent_matrix[j, i], 3:nc] * drop(exp(beta %*% pseudo_entries[index_recurrent_matrix[j, i], 3:nc]))
    }
    S1t_recurrent[, i] <- sum
  }
  return(S1t_recurrent / length(Y))
}


#' S1(t) for Recurrent Event Sub-Model with Cross-Fitting
#'
#' Computes the weighted first moment S1(t) at each recurrent event time
#' for the recurrence sub-model, using fold-specific beta coefficients
#' for cross-fitting.
#'
#' @param Y Vector of composite censoring/death times per subject.
#' @param wt_recurrent_subject Weight matrix at recurrent event times.
#' @param tr Vector of recurrent event times.
#' @param index_recurrent_matrix Matrix of pseudo entry indices for recurrent times.
#' @param pseudo_entries Sorted pseudo data set.
#' @param beta_mat Matrix of beta coefficients (rows = folds).
#' @param CV_map Two-column matrix mapping fold index to subject index.
#'
#' @return Matrix of S1(t) values (rows = covariates, cols = recurrent times).
#'
#' @keywords internal
jfm_s1t_recurrent_cf <- function(Y, wt_recurrent_subject, tr, index_recurrent_matrix, pseudo_entries, beta_mat, CV_map) {
  tr_ranked <- tr[order(tr)]
  nc <- ncol(pseudo_entries)
  # define a vector to store S1t value for different death time
  S1t_recurrent <- matrix(NA, nrow = nc - 2, ncol = length(tr))
  for (i in 1:length(tr)) {
    sum <- 0
    for (j in 1:length(Y)) {
      fold_j <- CV_map[CV_map[, 2] == j, 1]
      beta <- beta_mat[fold_j, ]
      sum <- sum + jfm_yt(Y[j], tr_ranked[i]) * wt_recurrent_subject[j, i] * pseudo_entries[index_recurrent_matrix[j, i], 3:nc] * drop(exp(beta %*% pseudo_entries[index_recurrent_matrix[j, i], 3:nc]))
    }
    S1t_recurrent[, i] <- sum
  }
  return(S1t_recurrent / length(Y))
}


#' S0(t) for Death Sub-Model
#'
#' Computes the weighted zeroth moment S0(t) at each death event time
#' for the death sub-model.
#'
#' @param Y Vector of composite censoring/death times per subject.
#' @param wt_matrix Weight matrix from \code{jfm_wt_death}.
#' @param td Vector of death event times.
#' @param index_death_matrix Matrix of pseudo entry indices for death times.
#' @param pseudo_entries Sorted pseudo data set.
#' @param alpha Coefficient vector for the death sub-model.
#'
#' @return Numeric vector of S0(t) values, one per death time.
#'
#' @keywords internal
jfm_s0t_death <- function(Y, wt_matrix, td, index_death_matrix, pseudo_entries, alpha) {
  td_ranked <- td[order(td)]
  # define a vector to store S1t value for different death time
  S0t_death <- numeric(length(td))
  nc <- ncol(pseudo_entries)
  for (i in 1:length(td)) {
    sum <- 0
    for (j in 1:length(Y)) {
      sum <- sum + jfm_yt(Y[j], td_ranked[i]) * wt_matrix[j, i] * drop(exp(alpha %*% pseudo_entries[index_death_matrix[j, i], 3:nc]))
    }
    S0t_death[i] <- sum
  }
  return(S0t_death / length(Y))
}


#' S0(t) for Death Sub-Model with Cross-Fitting
#'
#' Computes the weighted zeroth moment S0(t) at each death event time
#' for the death sub-model, using fold-specific alpha coefficients
#' for cross-fitting.
#'
#' @param Y Vector of composite censoring/death times per subject.
#' @param wt_matrix Weight matrix from \code{jfm_wt_death}.
#' @param td Vector of death event times.
#' @param index_death_matrix Matrix of pseudo entry indices for death times.
#' @param pseudo_entries Sorted pseudo data set.
#' @param alpha_mat Matrix of alpha coefficients (rows = folds).
#' @param CV_map Two-column matrix mapping fold index to subject index.
#'
#' @return Numeric vector of S0(t) values, one per death time.
#'
#' @keywords internal
jfm_s0t_death_cf <- function(Y, wt_matrix, td, index_death_matrix, pseudo_entries, alpha_mat, CV_map) {
  td_ranked <- td[order(td)]
  # define a vector to store S1t value for different death time
  S0t_death <- numeric(length(td))
  nc <- ncol(pseudo_entries)
  for (i in 1:length(td)) {
    sum <- 0
    for (j in 1:length(Y)) {
      fold_j <- CV_map[CV_map[, 2] == j, 1]
      alpha <- alpha_mat[fold_j, ]
      sum <- sum + jfm_yt(Y[j], td_ranked[i]) * wt_matrix[j, i] * drop(exp(alpha %*% pseudo_entries[index_death_matrix[j, i], 3:nc]))
    }
    S0t_death[i] <- sum
  }
  return(S0t_death / length(Y))
}


#' S0(t) for Recurrent Event Sub-Model
#'
#' Computes the weighted zeroth moment S0(t) at each recurrent event time
#' for the recurrence sub-model.
#'
#' @param Y Vector of composite censoring/death times per subject.
#' @param wt_recurrent_subject Weight matrix at recurrent event times.
#' @param tr Vector of recurrent event times.
#' @param index_recurrent_matrix Matrix of pseudo entry indices for recurrent times.
#' @param pseudo_entries Sorted pseudo data set.
#' @param beta Coefficient vector for the recurrent event sub-model.
#'
#' @return Numeric vector of S0(t) values, one per recurrent event time.
#'
#' @keywords internal
jfm_s0t_recurrent <- function(Y, wt_recurrent_subject, tr, index_recurrent_matrix, pseudo_entries, beta) {
  tr_ranked <- tr[order(tr)]
  # define a vector to store S1t value for different death time
  S0t_recurrent <- numeric(length(tr))
  nc <- ncol(pseudo_entries)
  for (i in 1:length(tr)) {
    sum <- 0
    for (j in 1:length(Y)) {
      sum <- sum + jfm_yt(Y[j], tr_ranked[i]) * wt_recurrent_subject[j, i] * drop(exp(beta %*% pseudo_entries[index_recurrent_matrix[j, i], 3:nc]))
    }
    S0t_recurrent[i] <- sum
  }
  return(S0t_recurrent / length(Y))
}


#' S0(t) for Recurrent Event Sub-Model with Cross-Fitting
#'
#' Computes the weighted zeroth moment S0(t) at each recurrent event time
#' for the recurrence sub-model, using fold-specific beta coefficients
#' for cross-fitting.
#'
#' @param Y Vector of composite censoring/death times per subject.
#' @param wt_recurrent_subject Weight matrix at recurrent event times.
#' @param tr Vector of recurrent event times.
#' @param index_recurrent_matrix Matrix of pseudo entry indices for recurrent times.
#' @param pseudo_entries Sorted pseudo data set.
#' @param beta_mat Matrix of beta coefficients (rows = folds).
#' @param CV_map Two-column matrix mapping fold index to subject index.
#'
#' @return Numeric vector of S0(t) values, one per recurrent event time.
#'
#' @keywords internal
jfm_s0t_recurrent_cf <- function(Y, wt_recurrent_subject, tr, index_recurrent_matrix, pseudo_entries, beta_mat, CV_map) {
  tr_ranked <- tr[order(tr)]
  # define a vector to store S1t value for different death time
  S0t_recurrent <- numeric(length(tr))
  nc <- ncol(pseudo_entries)
  for (i in 1:length(tr)) {
    sum <- 0
    for (j in 1:length(Y)) {
      fold_j <- CV_map[CV_map[, 2] == j, 1]
      beta <- beta_mat[fold_j, ]
      sum <- sum + jfm_yt(Y[j], tr_ranked[i]) * wt_recurrent_subject[j, i] * drop(exp(beta %*% pseudo_entries[index_recurrent_matrix[j, i], 3:nc]))
    }
    S0t_recurrent[i] <- sum
  }
  return(S0t_recurrent / length(Y))
}


#' Score Equation U1 for Beta (Recurrence)
#'
#' Evaluates the score equation for the recurrence coefficient vector beta.
#'
#' @param pseudo_entries Sorted pseudo data set.
#' @param index_recurrent_matrix Matrix of pseudo entry indices for recurrent times.
#' @param tr_id Reordered subject IDs for recurrent events.
#' @param S1t_re S1(t) matrix for recurrence.
#' @param S0t_re S0(t) vector for recurrence.
#'
#' @return Numeric vector of score values, one per covariate.
#'
#' @keywords internal
jfm_score_beta <- function(pseudo_entries, index_recurrent_matrix, tr_id, S1t_re, S0t_re) {
  U1 <- 0
  nc <- ncol(pseudo_entries)
  for (i in 1:length(S0t_re)) {
    U1 <- U1 + pseudo_entries[index_recurrent_matrix[tr_id[i], i], 3:nc] - S1t_re[, i] / S0t_re[i]
  }
  return(U1)
}


#' Score Equation U2 for Alpha (Death)
#'
#' Evaluates the score equation for the death coefficient vector alpha.
#'
#' @param pseudo_entries Sorted pseudo data set.
#' @param index_death_matrix Matrix of pseudo entry indices for death times.
#' @param td_id Reordered subject IDs for death events.
#' @param S1t_de S1(t) matrix for death.
#' @param S0t_de S0(t) vector for death.
#'
#' @return Numeric vector of score values, one per covariate.
#'
#' @keywords internal
jfm_score_alpha <- function(pseudo_entries, index_death_matrix, td_id, S1t_de, S0t_de) {
  U2 <- 0
  nc <- ncol(pseudo_entries)
  for (i in 1:length(S0t_de)) {
    U2 <- U2 + pseudo_entries[index_death_matrix[td_id[i], i], 3:nc] - S1t_de[, i] / S0t_de[i]
  }
  return(U2)
}


#' Estimating Equation U4 for Baseline Death Hazard
#'
#' Evaluates the estimating equation for the baseline hazard of the
#' death sub-model at each death event time.
#'
#' @param td Vector of death event times.
#' @param d_td Table of death time frequencies.
#' @param n Number of subjects.
#' @param S0t_de S0(t) vector for death.
#' @param lambda0_d Baseline hazard point masses for death.
#'
#' @return Named numeric vector of U4 values, one per death time.
#'
#' @keywords internal
jfm_u4 <- function(td, d_td, n, S0t_de, lambda0_d) {
  map_td <- order(td)
  d_td_ranked <- d_td[map_td]
  lambda0_d_ranked <- lambda0_d # lambda0_d_ranked=lambda0_d[map_td]
  U4 <- d_td_ranked - n * S0t_de * lambda0_d_ranked
  names(U4) <- 1:length(td)
  return(U4)
}


#' Solve for Baseline Death Hazard lambda0_d
#'
#' Computes the closed-form solution for the baseline hazard point masses
#' of the death sub-model.
#'
#' @param td Vector of death event times.
#' @param d_td Table of death time frequencies.
#' @param n Number of subjects.
#' @param S0t_de S0(t) vector for death.
#'
#' @return Named numeric vector of baseline hazard point masses for death.
#'
#' @keywords internal
jfm_lambda0d_solution <- function(td, d_td, n, S0t_de) {
  map_td <- order(td)
  d_td_ranked <- d_td[map_td]
  lambda0_d_ranked <- d_td_ranked / (n * S0t_de)
  names(lambda0_d_ranked) <- 1:length(td)
  return(lambda0_d_ranked)
}


#' Solve for Baseline Recurrence Hazard lambda0_r
#'
#' Computes the closed-form solution for the baseline hazard point masses
#' of the recurrence sub-model.
#'
#' @param tr Vector of recurrent event times.
#' @param d_tr Table of recurrent event frequencies.
#' @param n Number of subjects.
#' @param S0t_re S0(t) vector for recurrence.
#'
#' @return Named numeric vector of baseline hazard point masses for recurrence.
#'
#' @keywords internal
jfm_lambda0r_solution <- function(tr, d_tr, n, S0t_re) {
  map_tr <- order(tr)
  d_tr_ranked <- d_tr[map_tr]
  lambda0_r_ranked <- d_tr_ranked / (n * S0t_re)
  names(lambda0_r_ranked) <- 1:length(tr)
  return(lambda0_r_ranked)
}


#' Estimating Equation U5 for Baseline Recurrence Hazard
#'
#' Evaluates the estimating equation for the baseline hazard of the
#' recurrence sub-model at each recurrent event time.
#'
#' @param tr Vector of recurrent event times.
#' @param d_tr Table of recurrent event frequencies.
#' @param n Number of subjects.
#' @param S0t_re S0(t) vector for recurrence.
#' @param lambda0_r Baseline hazard point masses for recurrence.
#'
#' @return Named numeric vector of U5 values, one per recurrent event time.
#'
#' @keywords internal
jfm_u5 <- function(tr, d_tr, n, S0t_re, lambda0_r) {
  map_tr <- order(tr)
  d_tr_ranked <- d_tr[map_tr]
  lambda0_r_ranked <- lambda0_r # lambda0_r_ranked=lambda0_r[map_tr]
  U5 <- d_tr_ranked - n * S0t_re * lambda0_r_ranked
  names(U5) <- 1:length(tr)
  return(U5)
}


#' Theta Estimating Equation for the Stagewise Algorithm
#'
#' Self-contained version of the theta estimating equation that rebuilds
#' all intermediate quantities (wt_matrix, r2i, G-bar) internally, used
#' as the objective in a root-finding procedure for theta.
#'
#' @param theta Frailty variance parameter.
#' @param t.start Vector of interval start times.
#' @param I Vector of subject indicators for each pseudo entry.
#' @param Z List of covariate matrices, one per subject.
#' @param alpha Coefficient vector for the death sub-model.
#' @param beta Coefficient vector for the recurrent event sub-model.
#' @param lambda0_r Baseline hazard point masses for recurrence.
#' @param lambda0_d Baseline hazard point masses for death.
#' @param td Vector of death event times.
#' @param tr Vector of recurrent event times.
#' @param tr.id Subject IDs corresponding to each recurrent event.
#' @param td.id Subject IDs corresponding to each death event.
#' @param Y Vector of composite censoring/death times per subject.
#' @param STATUS Vector of death indicators per subject.
#' @param list_recur List of recurrent event times per subject.
#' @param num_recur Integer vector of recurrent event counts per subject.
#'
#' @return Scalar value of the scaled estimating equation for theta.
#'
#' @keywords internal
jfm_est_theta_new <- function(theta, t.start, I, Z, alpha, beta, lambda0_r, lambda0_d, td, tr, tr.id, td.id, Y, STATUS, list_recur, num_recur) {
  ### Generate wt_matrix
  # construct the pseudo data set (Lj,Ij,zj)
  Z1 <- numeric(0)
  for (i in 1:length(Z)) {
    Z1 <- rbind(Z1, Z[[i]])
  }
  pseudo <- cbind(t.start, I, Z1)
  pseudo <- pseudo[order(pseudo[, 1]), ]
  # rank the death events and also reorder the corresponding point mass (retain the rank-index mapping)
  map_td <- order(td) # ranking to index of td
  td_ranked <- td[map_td]
  td.id <- td.id[map_td] # index of td to id

  lambda0_d_ranked <- lambda0_d # lambda0_d_ranked=lambda0_d[map_td]
  # define a vector to store the index of the pseudo entries whose L is the largest L <= the death time
  lar_L <- numeric(length(td))
  j <- 1
  L <- pseudo[1, 1]
  for (i in 1:length(td)) {
    while (L <= td_ranked[i] && j < nrow(pseudo)) {
      j <- j + 1
      L <- pseudo[j, 1]
    }
    if (L <= td_ranked[i] && j == nrow(pseudo)) {
      lar_L[i] <- j
    } else {
      j <- j - 1
      lar_L[i] <- j
    }
    L <- pseudo[j, 1]
  }
  # define a matrix to store the wt for different subject i and death time (column for death time and row for subject)
  wt_matrix <- matrix(NA, nrow = length(Z), ncol = length(td))
  # also define a matrix to store the index of the pseudo entry whose L is the largest L before each death time for each subject (column for death time and row for subject)
  index_death_matrix <- matrix(NA, nrow = length(Z), ncol = length(td))
  for (i in 1:length(Z)) {
    integral <- 0
    for (j in 1:length(td)) {
      k <- lar_L[j]
      while (pseudo[k, 2] != i) {
        k <- k - 1
      }
      index_death_matrix[i, j] <- k
      integral <- integral + lambda0_d_ranked[j] * exp(alpha %*% pseudo[k, 3:ncol(pseudo)])
      wij <- 1 / (1 + theta * integral)
      wt_matrix[i, j] <- wij
    }
  }
  td_id <- td.id

  ### generate the integral in r2i

  # rank the recurrent events and also reorder the corresponding point mass (retain the rank-index mapping)
  map_tr <- order(tr)
  tr_ranked <- tr[map_tr]
  tr_id <- tr.id[map_tr] # rank to id mapping
  lambda0_r_ranked <- lambda0_r # lambda0_r_ranked=lambda0_r[map_tr]
  # define a matrix to store the index of the pseudo entries whose L is the largest L <= the recurrent time for subject i
  lar_L <- matrix(NA, nrow = length(Z), ncol = length(tr))
  j <- 1
  L <- pseudo[1, 1]

  # Create empty vectors to store indices of pseudo entries which belongs to subject j and are before tr_ranked[i]
  vectors <- vector("list", length(Z))
  for (i in 1:length(tr)) {
    while (L < tr_ranked[i]) {
      vectors[[pseudo[j, 2]]] <- c(vectors[[pseudo[j, 2]]], j)
      j <- j + 1
      L <- pseudo[j, 1]
    }


    for (k in 1:length(Z)) {
      lar_L[k, i] <- vectors[[k]][length(vectors[[k]])]
    }
  }
  # define a matrix to store the integral in r2i for different subject i and recurrent time (column for recurrent time and row for subject)
  integral_matrix <- matrix(NA, nrow = length(Z), ncol = length(tr))
  for (i in 1:length(Z)) {
    integral <- 0
    for (j in 1:length(tr)) {
      k <- lar_L[i, j]
      integral <- integral + lambda0_r_ranked[j] * exp(beta %*% pseudo[k, 3:ncol(pseudo)])
      integral_matrix[i, j] <- integral
    }
  }

  index_recurrent_matrix <- lar_L

  ### generate r2i evalated at each death time

  # td_ranked=td[order(td)]
  tr_ranked <- tr[order(tr)]
  r2i_death_subject_matrix <- matrix(NA, nrow = nrow(integral_matrix), ncol = length(td))
  recu_index <- numeric(length(td))
  flag <- numeric(length(td))
  j <- 1
  for (k in 1:length(td)) {
    flag[k] <- 1
    death <- td_ranked[k]
    T <- tr_ranked[j]
    while (death > T && j < length(tr_ranked)) {
      j <- j + 1
      T <- tr_ranked[j]
    }
    if (death > T && j == length(tr_ranked)) {
      recu_index[k] <- j
    } else {
      if (death == T) {
        recu_index[k] <- j
        j <- j + 1
      } else {
        if (j == 1) {
          recu_index[k] <- j
          flag[k] <- 0
        } else {
          recu_index[k] <- j - 1
        }
      }
    }
  }

  for (i in 1:nrow(integral_matrix)) {
    for (j in 1:length(td)) {
      r2i_death_subject_matrix[i, j] <- flag[j] * wt_matrix[i, j] * integral_matrix[i, recu_index[j]]
    }
  }

  ### Generate Gbar evaluated at each death time
  Gt_bar_death_matrix <- jfm_gt_bar_death(r2i_death_subject_matrix, td, Y, STATUS, list_recur)
  est <- 0
  for (i in 1:length(td)) {
    est <- est + num_recur[td_id[i]] - (theta + 1) * Gt_bar_death_matrix[i] * r2i_death_subject_matrix[td_id[i], i]
  }
  return(est / length(Y))
}


#' Generate a Single Covariate Vector
#'
#' Generates n observations from a specified distribution for simulation
#' purposes.
#'
#' @param n Number of observations.
#' @param v.type Distribution type: 1 = normal, 2 = Bernoulli, 3 = Poisson.
#' @param v.coeffi Numeric vector of distribution parameters. For normal:
#'   c(mean, sd). For Bernoulli: c(prob). For Poisson: c(lambda).
#'
#' @return Numeric vector of length n.
#'
#' @keywords internal
jfm_v_gene <- function(n, v.type, v.coeffi) {
  ## Generate single covariate according to the different distributions and parameters
  if (v.type == 1) {
    res <- rnorm(n, mean = v.coeffi[1], sd = v.coeffi[2])
  }
  if (v.type == 2) {
    res <- ifelse(runif(n) < v.coeffi[1], 1, 0)
  }
  if (v.type == 3) {
    res <- rpois(n, lambda = v.coeffi[1])
  }
  res
}


#' Inverse CDF of Exponential Distribution
#'
#' Computes the quantile function (inverse CDF) of an exponential
#' distribution for simulation purposes.
#'
#' @param y Probability value(s) in (0, 1).
#' @param rate Rate parameter of the exponential distribution.
#'
#' @return Numeric value(s) corresponding to the inverse CDF.
#'
#' @keywords internal
jfm_exp_icdf <- function(y, rate) {
  return(-log(1 - y) / rate)
}


#' Estimating Equation U4 (New) for Numerically Solving Baseline Death Hazard
#'
#' Self-contained version of U4 that rebuilds the weight matrix and S0(t)
#' internally, used for numerically solving for lambda0_d when a
#' closed-form solution is not applicable.
#'
#' @param lambda0_d Baseline hazard point masses for death.
#' @param d_td Table of death time frequencies.
#' @param n Number of subjects.
#' @param Y Vector of composite censoring/death times per subject.
#' @param td Vector of death event times.
#' @param theta Frailty variance parameter.
#' @param alpha Coefficient vector for the death sub-model.
#' @param t.start Vector of interval start times.
#' @param I Vector of subject indicators for each pseudo entry.
#' @param Z List of covariate matrices, one per subject.
#' @param td.id Subject IDs corresponding to each death event.
#'
#' @return Named numeric vector of U4 values, one per death time.
#'
#' @keywords internal
jfm_u4_new <- function(lambda0_d, d_td, n, Y, td, theta, alpha, t.start, I, Z, td.id) {
  # generate S0t evaluated at each death time
  # construct wt matrix evaluated at each death time for each subject
  # construct the pseudo data set (Lj,Ij,zj)
  Z1 <- numeric(0)
  for (i in 1:length(Z)) {
    Z1 <- rbind(Z1, Z[[i]])
  }
  pseudo <- cbind(t.start, I, Z1)
  pseudo <- pseudo[order(pseudo[, 1]), ]
  # rank the death events and also reorder the corresponding point mass (retain the rank-index mapping)
  map_td <- order(td) # ranking to index of td
  td_ranked <- td[map_td]
  td.id <- td.id[map_td] # index of td to id

  lambda0_d_ranked <- lambda0_d # lambda0_d_ranked=lambda0_d[map_td]
  # define a vector to store the index of the pseudo entries whose L is the largest L <= the death time
  lar_L <- numeric(length(td))
  j <- 1
  L <- pseudo[1, 1]
  for (i in 1:length(td)) {
    while (L <= td_ranked[i] && j < nrow(pseudo)) {
      j <- j + 1
      L <- pseudo[j, 1]
    }
    if (L <= td_ranked[i] && j == nrow(pseudo)) {
      lar_L[i] <- j
    } else {
      j <- j - 1
      lar_L[i] <- j
    }
    L <- pseudo[j, 1]
  }
  # define a matrix to store the wt for different subject i and death time (column for death time and row for subject)
  wt_matrix <- matrix(NA, nrow = length(Z), ncol = length(td))
  # also define a matrix to store the index of the pseudo entry whose L is the largest L before each death time for each subject (column for death time and row for subject)
  index_death_matrix <- matrix(NA, nrow = length(Z), ncol = length(td))
  for (i in 1:length(Z)) {
    integral <- 0
    for (j in 1:length(td)) {
      k <- lar_L[j]
      while (pseudo[k, 2] != i) {
        k <- k - 1
      }
      index_death_matrix[i, j] <- k
      integral <- integral + lambda0_d_ranked[j] * exp(alpha %*% pseudo[k, 3:ncol(pseudo)])
      wij <- 1 / (1 + theta * integral)
      wt_matrix[i, j] <- wij
    }
  }

  td_ranked <- td[order(td)]
  # define a vector to store S1t value for different death time
  S0t_death <- numeric(length(td))
  nc <- ncol(pseudo)
  for (i in 1:length(td)) {
    sum <- 0
    for (j in 1:length(Y)) {
      sum <- sum + jfm_yt(Y[j], td_ranked[i]) * wt_matrix[j, i] * drop(exp(alpha %*% pseudo[index_death_matrix[j, i], 3:nc]))
    }
    S0t_death[i] <- sum
  }
  S0t_de <- S0t_death / length(Y)


  map_td <- order(td)
  d_td_ranked <- d_td[map_td]
  lambda0_d_ranked <- lambda0_d # lambda0_d_ranked=lambda0_d[map_td]
  U4 <- d_td_ranked - n * S0t_de * lambda0_d_ranked
  names(U4) <- 1:length(td)
  return(U4)
}
