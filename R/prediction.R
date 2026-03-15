# ============================================================================
# Prediction helpers for time-dependent AUC and CIF markers (JFM + JSCM)
# ============================================================================


# Internal: Nelson-Aalen baseline hazard for JSCM on the accelerated time scale.
#
# For the JSCM the hazard is λ(t|Z_i) = e^{α^T Z_i} λ_0(t · e^{α^T Z_i}),
# so the cumulative hazard for subject i at time t is
#   Λ(t|Z_i) = Λ_0(t · e^{α^T Z_i})
# and S(t|Z_i) = exp(-Λ_0(t · e^{α^T Z_i})).
#
# The Nelson-Aalen estimator is computed by pooling all (accelerated) event
# times and using each subject's accelerated follow-up time as the censoring
# boundary for the risk set.
.jscm_baseline <- function(data, alpha, beta) {
  p        <- sum(startsWith(names(data), "x"))
  cov_cols <- paste0("x", seq_len(p))

  # Per-subject row (first occurrence) — covariates are time-invariant
  sub_rows <- !duplicated(data$id)
  ids      <- data$id[sub_rows]
  Z        <- as.matrix(data[sub_rows, cov_cols, drop = FALSE])
  id_idx   <- match(data$id, ids)   # row index of each data row in Z

  lp_re    <- as.numeric(Z %*% alpha)
  lp_de    <- as.numeric(Z %*% beta)
  accel_re <- exp(lp_re)
  accel_de <- exp(lp_de)

  # Follow-up time per subject: max t.stop
  fu_raw   <- tapply(data$t.stop, data$id, max)[as.character(ids)]
  fu_re    <- fu_raw * accel_re   # on accelerated recurrence scale
  fu_de    <- fu_raw * accel_de   # on accelerated death scale

  # --- Recurrence (all events with event == 1) ---
  rec      <- data[data$event == 1, ]
  t_r_acc  <- rec$t.stop * accel_re[id_idx[data$event == 1]]
  ord_r    <- order(t_r_acc)
  t_r_s    <- t_r_acc[ord_r]
  # Risk set at each (sorted) accelerated event time
  denom_r  <- vapply(t_r_s, function(t) sum(fu_re >= t), integer(1L))
  denom_r  <- pmax(denom_r, 1L)

  # --- Death (last observation with status == 1) ---
  last_rows <- !duplicated(data$id, fromLast = TRUE)
  t_de_all  <- data$t.stop[last_rows]
  st_de     <- data$status[last_rows]
  idx_de    <- which(last_rows)
  t_d_acc   <- t_de_all * accel_de[id_idx[last_rows]]

  death_rows <- st_de == 1
  t_d_s_all  <- t_d_acc[death_rows]
  ord_d      <- order(t_d_s_all)
  t_d_s      <- t_d_s_all[ord_d]
  denom_d    <- vapply(t_d_s, function(t) sum(fu_de >= t), integer(1L))
  denom_d    <- pmax(denom_d, 1L)

  list(
    t_r      = t_r_s,
    lambda0_r = 1 / denom_r,
    t_d      = t_d_s,
    lambda0_d = 1 / denom_d
  )
}


#' Cumulative Baseline Hazard Step Functions
#'
#' Evaluates the cumulative baseline hazards for readmission and death at
#' arbitrary time points.  For JFM, Breslow-type estimators are used.  For
#' JSCM, Nelson-Aalen estimators on the accelerated time scale are used.
#'
#' @param object An object of class \code{"swjm_cv"}.
#' @param times Numeric vector of evaluation times. If \code{NULL}, the
#'   observed event times stored in the fit are used.
#' @param which Character. Which baseline hazard(s) to return:
#'   \code{"both"} (default), \code{"readmission"}, or \code{"death"}.
#'
#' @return A data frame with column \code{time} and, depending on
#'   \code{which}, \code{cumhaz_readmission} and/or \code{cumhaz_death}.
#'
#' @examples
#' \donttest{
#' dat <- generate_data(n = 50, p = 5, scenario = 1, model = "jfm")
#' cv  <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
#'                     max_iter = 100)
#' bh  <- baseline_hazard(cv)
#' head(bh)
#' }
#'
#' @export
baseline_hazard <- function(object,
                             times = NULL,
                             which = c("both", "readmission", "death")) {
  which <- match.arg(which)
  if (!inherits(object, "swjm_cv"))
    stop("'object' must be of class 'swjm_cv'")
  if (!object$model %in% c("jfm", "jscm"))
    stop("baseline_hazard() requires model = 'jfm' or 'jscm'")
  if (is.null(object$baseline))
    stop("No baseline hazard stored in this fit")

  bl        <- object$baseline
  # JFM stores event times as tr/td; JSCM (accelerated scale) as t_r/t_d
  tr_sorted <- if (object$model == "jfm") bl$tr else bl$t_r
  td_sorted <- if (object$model == "jfm") bl$td else bl$t_d

  if (is.null(times))
    times <- sort(unique(c(tr_sorted, td_sorted)))

  # Prepend 0 so index 0 (before first event) maps to 0
  cum_r    <- c(0, cumsum(bl$lambda0_r))
  cum_d    <- c(0, cumsum(bl$lambda0_d))
  idx_r    <- findInterval(times, tr_sorted) + 1L
  idx_d    <- findInterval(times, td_sorted) + 1L
  cumhaz_r <- cum_r[idx_r]
  cumhaz_d <- cum_d[idx_d]

  out <- data.frame(time = times)
  if (which %in% c("both", "readmission")) out$cumhaz_readmission <- cumhaz_r
  if (which %in% c("both", "death"))       out$cumhaz_death       <- cumhaz_d
  out
}


#' Cumulative Hazard from Estimated Baseline Hazards (JFM)
#'
#' Computes subject-specific cumulative hazards for recurrent events using
#' estimated baseline hazard point masses and the pseudo data approach.
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
#'   \item{integral_matrix}{Matrix of cumulative hazard values
#'     (rows = subjects, cols = recurrent times).}
#'   \item{index_recurrent_matrix}{Matrix of pseudo entry indices.}
#'   \item{tr_id}{Reordered subject IDs for recurrent events.}
#'
#' @keywords internal
cumulative_hazard_jfm <- function(t.start, I, Z, beta, tr, lambda0_r, tr.id) {
  Z1 <- numeric(0)
  for (i in seq_along(Z)) {
    Z1 <- rbind(Z1, Z[[i]])
  }
  pseudo <- cbind(t.start, I, Z1)
  pseudo <- pseudo[order(pseudo[, 1]), ]

  map_tr <- order(tr)
  tr_ranked <- tr[map_tr]
  tr_id <- tr.id[map_tr]
  lambda0_r_ranked <- lambda0_r

  lar_L <- matrix(NA, nrow = length(Z), ncol = length(tr))
  j <- 1L
  L <- pseudo[1, 1]

  vectors <- vector("list", length(Z))
  for (i in seq_along(tr)) {
    while (L < tr_ranked[i] && j < nrow(pseudo)) {
      vectors[[pseudo[j, 2]]] <- c(vectors[[pseudo[j, 2]]], j)
      j <- j + 1L
      L <- pseudo[j, 1]
    }
    if (L < tr_ranked[i] && j == nrow(pseudo)) {
      vectors[[pseudo[j, 2]]] <- c(vectors[[pseudo[j, 2]]], j)
    }
    for (k in seq_along(Z)) {
      lar_L[k, i] <- vectors[[k]][length(vectors[[k]])]
    }
  }

  integral_matrix <- matrix(NA, nrow = length(Z), ncol = length(tr))
  for (i in seq_along(Z)) {
    integral <- 0
    for (j2 in seq_along(tr)) {
      k <- lar_L[i, j2]
      integral <- integral + lambda0_r_ranked[j2] *
        exp(beta %*% pseudo[k, 3:ncol(pseudo)])
      integral_matrix[i, j2] <- integral
    }
  }
  list(integral_matrix = integral_matrix,
       index_recurrent_matrix = lar_L,
       tr_id = tr_id)
}


#' Cumulative Hazard Under True Baseline Hazard (JFM)
#'
#' Computes the cumulative hazard at time \code{x} under piecewise-constant
#' baseline hazard with known cut points and covariate values.
#'
#' @param x Scalar time at which to evaluate.
#' @param cut_points Numeric vector of cut point times.
#' @param z_values Matrix of covariate values per interval.
#' @param beta Coefficient vector.
#'
#' @return Scalar cumulative hazard value.
#'
#' @keywords internal
cumulative_hazard_true_jfm <- function(x, cut_points, z_values, beta) {
  integral <- 0
  for (i in seq_len(nrow(z_values))) {
    if (i == 1) {
      start <- 0
      end <- min(x, cut_points[1])
    } else if (i <= length(cut_points)) {
      start <- cut_points[i - 1]
      end <- min(x, cut_points[i])
    } else {
      start <- cut_points[length(cut_points)]
      end <- x
    }
    if (start >= x) break
    integral <- integral + as.numeric(exp(beta %*% z_values[i, ])) *
      (end - start)
  }
  integral
}


#' Cause-Specific CIF Markers (JFM)
#'
#' Computes cause-specific cumulative incidence function (CIF) markers
#' for recurrent events under the joint frailty model, using estimated
#' baseline hazards.
#'
#' @param alpha Coefficient vector for the recurrence sub-model
#'   (first p elements of theta).
#' @param beta Coefficient vector for the death sub-model
#'   (second p elements of theta).
#' @param tr Vector of recurrent event times.
#' @param td Vector of death event times.
#' @param lambda0_r Baseline hazard point masses for recurrence.
#' @param lambda0_d Baseline hazard point masses for death.
#' @param Z_base Matrix of baseline covariate values (one row per subject).
#'
#' @return A list with components:
#'   \item{cif_marker}{Matrix of CIF values (rows = subjects, cols = event times).}
#'   \item{linear_marker}{Numeric vector of linear predictors for recurrence.}
#'
#' @keywords internal
markers_cause_specific <- function(alpha, beta, tr, td, lambda0_r, lambda0_d,
                                   Z_base) {
  # Cumulative baseline hazard of recurrence
  cum_lambda0_r <- cumsum(lambda0_r)

  # Cumulative baseline hazard of death
  cum_lambda0_d <- cumsum(lambda0_d)

  # Map death cumulative hazard to recurrence event times
  map_tr <- order(tr)
  tr_ranked <- tr[map_tr]
  map_td <- order(td)
  td_ranked <- td[map_td]

  idx <- findInterval(tr_ranked, td_ranked)
  cum_lambda0_d_r <- ifelse(idx == 0, 0, cum_lambda0_d[idx])

  # Compute per-subject CIF
  # alpha = recurrence coefficient, beta = death coefficient
  S_mat <- matrix(NA, nrow = nrow(Z_base), ncol = length(tr))
  cif_mat <- matrix(NA, nrow = nrow(Z_base), ncol = length(tr))
  eta_re_vec <- numeric(nrow(Z_base))

  for (i in seq_len(nrow(Z_base))) {
    eta_de <- drop(beta  %*% Z_base[i, ])   # death linear predictor
    eta_re <- drop(alpha %*% Z_base[i, ])   # recurrence linear predictor
    eta_re_vec[i] <- eta_re

    S_mat[i, ] <- exp(-(exp(eta_de) * cum_lambda0_d_r +
                           exp(eta_re) * cum_lambda0_r))
    cif_mat[i, ] <- exp(eta_re) * cumsum(S_mat[i, ] * lambda0_r)
  }

  list(cif_marker = cif_mat, linear_marker = eta_re_vec)
}


#' Cause-Specific CIF Markers Under True Parameters (JFM)
#'
#' Computes the CIF at a single time point under constant baseline hazards.
#'
#' @param t Scalar time at which to evaluate.
#' @param alpha Coefficient vector for the recurrence sub-model
#'   (first p elements of theta).
#' @param beta Coefficient vector for the death sub-model
#'   (second p elements of theta).
#' @param lambda0_r Scalar baseline hazard rate for recurrence.
#' @param lambda0_d Scalar baseline hazard rate for death.
#' @param Z_base Matrix of baseline covariate values (one row per subject).
#'
#' @return A list with components:
#'   \item{cif_marker}{Numeric vector of CIF values per subject.}
#'   \item{linear_marker}{Numeric vector of linear predictors for recurrence.}
#'
#' @keywords internal
markers_cause_specific_true <- function(t, alpha, beta, lambda0_r, lambda0_d,
                                        Z_base) {
  cif_marker <- numeric(nrow(Z_base))
  eta_re_vec <- numeric(nrow(Z_base))

  for (i in seq_len(nrow(Z_base))) {
    eta_de <- drop(beta  %*% Z_base[i, ])   # death linear predictor
    eta_re <- drop(alpha %*% Z_base[i, ])   # recurrence linear predictor
    eta_re_vec[i] <- eta_re

    A <- lambda0_r * exp(eta_re) + lambda0_d * exp(eta_de)
    cif_marker[i] <- (exp(eta_re) / A) * (1 - exp(-A * t))
  }

  list(cif_marker = cif_marker, linear_marker = eta_re_vec)
}


#' Generate Baseline Hazard for Recurrence (JFM)
#'
#' Computes the Breslow-type baseline hazard for recurrence from
#' a fitted alpha coefficient vector.
#'
#' @param Data2 A data frame in recurrent-event format.
#' @param alpha Coefficient vector for the recurrence sub-model
#'   (first p elements of theta).
#'
#' @return Named numeric vector of baseline hazard point masses.
#'
#' @keywords internal
lambda0_r_gen <- function(Data2, alpha) {
  dc <- extract_data_components(Data2)
  Z <- dc$Z; n <- dc$n; tr <- dc$tr; tr.id <- dc$tr.id; d_tr <- dc$d_tr
  td <- dc$td; td.id <- dc$td.id; Y <- dc$Y
  t.start <- dc$t.start; I <- dc$I
  p <- dc$p

  lambda0_d_init <- rep(0.2, length(td))
  theta <- 1
  beta0 <- numeric(p)

  A <- jfm_wt_death(theta, beta0, t.start, I, Z, td, lambda0_d_init, td.id)
  pseudo_entries <- A$pseudo_entries

  diff_tr1 <- diff(c(0, tr[order(tr)]))
  lambda0_r_init <- diff_tr1 * 1
  B <- jfm_r2i_integral(t.start, I, Z, alpha, tr, lambda0_r_init, tr.id)
  index_recurrent_matrix <- B$index_recurrent_matrix

  wt_recurrent_subject <- matrix(1, nrow = n, ncol = length(tr))
  S0t_re <- jfm_s0t_recurrent(Y, wt_recurrent_subject, tr,
                               index_recurrent_matrix, pseudo_entries, alpha)
  jfm_lambda0r_solution(tr, d_tr, n, S0t_re)
}


#' Generate Both Baseline Hazards (JFM)
#'
#' Computes Breslow-type baseline hazards for both death and recurrence.
#'
#' @param Data2 A data frame in recurrent-event format.
#' @param alpha Coefficient vector for the recurrence sub-model
#'   (first p elements of theta).
#' @param beta Coefficient vector for the death sub-model
#'   (second p elements of theta).
#'
#' @return A list with components \code{lambda0_r} and \code{lambda0_d}.
#'
#' @keywords internal
lambda0_gen <- function(Data2, alpha, beta) {
  dc <- extract_data_components(Data2)
  Z <- dc$Z; n <- dc$n; td <- dc$td; td.id <- dc$td.id; d_td <- dc$d_td
  tr <- dc$tr; tr.id <- dc$tr.id; d_tr <- dc$d_tr
  Y <- dc$Y; t.start <- dc$t.start; I <- dc$I
  p <- dc$p

  lambda0_d_init <- rep(0.2, length(td))
  theta <- 1
  beta0 <- numeric(p)

  A <- jfm_wt_death(theta, beta0, t.start, I, Z, td,
                     lambda0_d_init, td.id)
  td_id <- A$td_id
  index_death_matrix <- A$index_death_matrix
  pseudo_entries <- A$pseudo_entries
  wt_matrix <- matrix(1, nrow = n, ncol = length(td))

  # alpha = recurrence, beta = death
  S0t_de <- jfm_s0t_death(Y, wt_matrix, td, index_death_matrix,
                           pseudo_entries, beta)
  lambda0_d <- jfm_lambda0d_solution(td, d_td, n, S0t_de)

  diff_tr1 <- diff(c(0, tr[order(tr)]))
  lambda0_r_init <- diff_tr1 * 1
  B <- jfm_r2i_integral(t.start, I, Z, alpha, tr, lambda0_r_init, tr.id)
  index_recurrent_matrix <- B$index_recurrent_matrix

  wt_recurrent_subject <- matrix(1, nrow = n, ncol = length(tr))
  S0t_re <- jfm_s0t_recurrent(Y, wt_recurrent_subject, tr,
                               index_recurrent_matrix, pseudo_entries, alpha)
  lambda0_r <- jfm_lambda0r_solution(tr, d_tr, n, S0t_re)

  list(lambda0_r = lambda0_r, lambda0_d = lambda0_d)
}
