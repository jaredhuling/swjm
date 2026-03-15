# ============================================================================
# Stagewise fitting for both JFM and JSCM
# ============================================================================


#' Fit a Stagewise Regularization Path
#'
#' Unified interface for stagewise variable selection in joint models
#' of recurrent events and terminal events. Dispatches to model-specific
#' implementations for the Joint Frailty Model (JFM) or Joint Scale-Change
#' Model (JSCM).
#'
#' @param data A data frame in recurrent-event format with columns
#'   \code{id}, \code{t.start}, \code{t.stop}, \code{event}, \code{status},
#'   and covariate columns \code{x1}, \ldots, \code{xp}.
#' @param model Character. Either \code{"jfm"} or \code{"jscm"}.
#' @param penalty Character. One of \code{"coop"} (cooperative lasso),
#'   \code{"lasso"}, or \code{"group"} (group lasso).
#' @param eps Numeric. Step size for the stagewise update. If \code{NULL},
#'   uses adaptive step size.
#' @param max_iter Integer. Maximum number of stagewise iterations.
#' @param pp Integer. Early-stopping block size: algorithm checks every
#'   \code{pp} iterations if fewer than 3 unique coordinates were updated.
#'
#' @return An object of class \code{"swjm_path"}, a list with components:
#'   \describe{
#'     \item{k}{Number of iterations performed.}
#'     \item{theta}{Matrix of coefficient paths (\code{2p} rows by
#'       \code{k+1} columns).}
#'     \item{lambda}{Numeric vector of penalty parameter approximations
#'       at each iteration.}
#'     \item{norm}{Numeric vector of penalty norm values along the path.}
#'     \item{model}{Character, the model used.}
#'     \item{penalty}{Character, the penalty used.}
#'   }
#'
#' @examples
#' \donttest{
#' dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
#' fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
#'                      max_iter = 100)
#' fit
#' }
#'
#' @export
stagewise_fit <- function(data, model = c("jfm", "jscm"),
                          penalty = c("coop", "lasso", "group"),
                          eps = NULL, max_iter = NULL, pp = NULL) {
  model <- match.arg(model)
  penalty <- match.arg(penalty)

  p <- ncol(data) - 5L
  initial_alpha <- numeric(p)
  initial_beta <- numeric(p)

  if (model == "jfm") {
    if (is.null(eps)) eps <- 0.1
    if (is.null(max_iter)) max_iter <- 5000L
    if (is.null(pp)) pp <- max_iter   # disable early stopping by default for JFM
    result <- stagewise_jfm(initial_alpha, initial_beta, data, penalty,
                            eps1 = 1e-6, adap = 1L, eps2 = eps,
                            iter = max_iter, pp = pp)
  } else {
    if (is.null(eps)) eps <- 0.01
    if (is.null(max_iter)) max_iter <- 5000L
    if (is.null(pp)) pp <- max_iter   # disable early stopping by default for JSCM
    result <- stagewise_jscm(initial_alpha, initial_beta, data, penalty,
                             eps1 = 1e-6, adap = 1L, eps2 = eps,
                             iter = max_iter, pp = pp)
  }

  structure(
    list(
      k = result$k,
      theta = result$theta_update,
      alpha = result$theta_update[1:p, , drop = FALSE],
      beta  = result$theta_update[(p + 1):(2 * p), , drop = FALSE],
      lambda = result$lambda,
      norm = result$normK_update,
      norm_g1 = result$norm2_g1_update,
      norm_g2 = result$norm2_g2_update,
      g1 = result$g1_update,
      g2 = result$g2_update,
      model = model,
      penalty = penalty,
      p = p
    ),
    class = "swjm_path"
  )
}


#' @export
print.swjm_path <- function(x, ...) {
  cat("Stagewise path (", x$model, "/", x$penalty, ")\n\n", sep = "")
  cat(sprintf("  %-26s %d\n", "Covariates (p):", x$p))
  cat(sprintf("  %-26s %d\n", "Iterations:", x$k))
  cat(sprintf("  %-26s [%.4g, %.4g]\n", "Lambda range:",
              min(x$lambda), max(x$lambda)))
  alpha_end <- x$alpha[, ncol(x$alpha)]
  beta_end  <- x$beta[, ncol(x$beta)]
  na <- sum(alpha_end != 0)
  nb <- sum(beta_end != 0)
  cat(sprintf("  %-26s %d readmission, %d death\n",
              "Active at final step:", na, nb))
  if (na > 0)
    cat(sprintf("    Readmission (alpha): %s\n",
                paste(which(alpha_end != 0), collapse = ", ")))
  if (nb > 0)
    cat(sprintf("    Death (beta):        %s\n",
                paste(which(beta_end != 0), collapse = ", ")))
  invisible(x)
}


#' @export
summary.swjm_path <- function(object, ...) {
  p   <- object$p
  lam <- object$lambda
  dec_idx   <- extract_decreasing_indices(lam)
  alpha_dec <- object$alpha[, dec_idx, drop = FALSE]
  beta_dec  <- object$beta[, dec_idx, drop = FALSE]

  alpha_end <- object$alpha[, ncol(object$alpha)]
  beta_end  <- object$beta[, ncol(object$beta)]

  # First entry step along the decreasing path (for sort order only)
  alpha_entry_step <- vapply(seq_len(p), function(j) {
    idx <- which(alpha_dec[j, ] != 0); if (length(idx)) idx[1L] else NA_integer_
  }, integer(1))
  beta_entry_step <- vapply(seq_len(p), function(j) {
    idx <- which(beta_dec[j, ] != 0); if (length(idx)) idx[1L] else NA_integer_
  }, integer(1))

  any_active   <- which(!is.na(alpha_entry_step) | !is.na(beta_entry_step))
  never_active <- setdiff(seq_len(p), any_active)

  cat(sprintf("Stagewise path (%s/%s)\n\n", object$model, object$penalty))
  cat(sprintf("  p = %d  |  %d iterations  |  lambda: [%.4g, %.4g]\n",
              p, object$k, min(lam), max(lam)))
  cat(sprintf("  Decreasing path: %d steps\n", length(dec_idx)))

  if (length(any_active) == 0L) {
    cat("\n  No variables entered the path.\n")
    return(invisible(object))
  }

  # Sort variables by first entry step
  first_entry <- pmin(
    ifelse(is.na(alpha_entry_step), .Machine$integer.max, alpha_entry_step),
    ifelse(is.na(beta_entry_step),  .Machine$integer.max, beta_entry_step)
  )
  ord <- any_active[order(first_entry[any_active])]

  # Determine variable type for each selected variable
  .var_type <- function(j) {
    a <- alpha_end[j]; b <- beta_end[j]
    if (a != 0 && b != 0) {
      if (sign(a) == sign(b)) "shared (+)" else "shared (\u2013)"
    } else if (a != 0) {
      "readmission only"
    } else {
      "death only"
    }
  }

  cat("\n  Path-end coefficients (nonzero variables):\n\n")
  cat(sprintf("  %-10s  %-10s  %-10s  %s\n",
              "Variable", "alpha", "beta", "Type"))
  cat("  ----------  ----------  ----------  ----------------\n")
  for (j in ord) {
    a_str <- if (alpha_end[j] != 0) sprintf("%+.4f", alpha_end[j]) else "     \u2014"
    b_str <- if (beta_end[j]  != 0) sprintf("%+.4f", beta_end[j])  else "     \u2014"
    cat(sprintf("  %-10s  %-10s  %-10s  %s\n",
                paste0("x", j), a_str, b_str, .var_type(j)))
  }

  if (length(never_active) > 0L)
    cat(sprintf("\n  Inactive: %s\n",
                paste(paste0("x", never_active), collapse = ", ")))

  invisible(object)
}


# --------------------------------------------------------------------------
# JFM stagewise implementation
# --------------------------------------------------------------------------

#' @keywords internal
stagewise_jfm <- function(initial_alpha, initial_beta, Data2, penalty,
                          eps1, adap, eps2, iter, pp = 200L) {
  p <- length(initial_alpha)
  dc <- extract_data_components(Data2)
  Z <- dc$Z; n <- dc$n; td <- dc$td; td.id <- dc$td.id; d_td <- dc$d_td
  tr <- dc$tr; tr.id <- dc$tr.id; d_tr <- dc$d_tr
  Y <- dc$Y; STATUS <- dc$STATUS
  list_recur <- dc$list_recur; num_recur <- dc$num_recur
  t.start <- dc$t.start; I <- dc$I

  # Initialize alpha (recurrence) and beta (death) -- consistent with JSCM
  alpha <- initial_alpha
  beta <- initial_beta

  # Initialize baseline hazards
  lambda0_d <- rep(0.2, length(td))
  theta <- 1
  A <- jfm_wt_death(theta, beta, t.start, I, Z, td, lambda0_d, td.id)
  td_id    <- A$td_id
  index_death_matrix    <- A$index_death_matrix
  pseudo_entries        <- A$pseudo_entries
  wt_matrix             <- matrix(1, nrow = n, ncol = length(td))

  # Pre-extract covariate block and sorted event times (used every iteration)
  Z_pseudo  <- pseudo_entries[, 3:ncol(pseudo_entries), drop = FALSE]
  td_sorted <- sort(td)
  # tr_id: 0-based subject index per recurrent event (for C++ score)
  # td_id: same for death -- computed after wt_death call above
  td_id0    <- td_id - 1L   # 0-based

  res_de <- jfm_s0s1_cpp(Y, wt_matrix, td_sorted, index_death_matrix,
                          Z_pseudo, beta)
  S0t_de <- res_de$S0t
  lambda0_d <- jfm_lambda0d_solution(td, d_td, n, S0t_de)

  diff_tr1 <- diff(c(0, tr[order(tr)]))
  lambda0_r <- diff_tr1 * 1

  B <- jfm_r2i_integral(t.start, I, Z, alpha, tr, lambda0_r, tr.id)
  index_recurrent_matrix <- B$index_recurrent_matrix
  tr_id  <- B$tr_id
  tr_id0 <- tr_id - 1L   # 0-based
  tr_sorted <- sort(tr)
  wt_recurrent_subject <- matrix(1, nrow = n, ncol = length(tr))

  res_re <- jfm_s0s1_cpp(Y, wt_recurrent_subject, tr_sorted,
                          index_recurrent_matrix, Z_pseudo, alpha)
  S0t_re <- res_re$S0t
  lambda0_r <- jfm_lambda0r_solution(tr, d_tr, n, S0t_re)

  # Stagewise loop
  k <- 0L
  thetaK <- c(initial_alpha, initial_beta)
  normK <- 0

  # Compute initial gradients using combined C++ calls
  g1 <- (-1) * jfm_score_cpp(index_recurrent_matrix, tr_id0,
                               Z_pseudo, res_re$S1t, S0t_re) / n
  g2 <- (-1) * jfm_score_cpp(index_death_matrix, td_id0,
                               Z_pseudo, res_de$S1t, S0t_de) / n

  # Gradient scaling: scale death (g2) up by max|g1|/max|g2|
  if (penalty %in% c("coop", "group")) {
    cc <- max(abs(g1)) / max(abs(g2))
  } else {
    cc <- sqrt(sum(g1^2)) / sqrt(sum(g2^2))
  }
  g2_origin <- g2
  g2 <- cc * g2

  # Storage
  theta_update <- matrix(NA, nrow = 2 * p, ncol = iter + 1L)
  theta_update[, 1] <- thetaK
  normK_update <- rep(NA_real_, iter + 1L)
  normK_update[1] <- normK
  norm2_g1_update <- rep(NA_real_, iter + 1L)
  norm2_g1_update[1] <- sqrt(sum(g1^2))
  norm2_g2_update <- rep(NA_real_, iter + 1L)
  norm2_g2_update[1] <- sqrt(sum(g2^2))
  lambda_vec <- rep(NA_real_, iter + 1L)
  g1_update <- matrix(NA, nrow = p, ncol = iter + 1L)
  g2_update <- matrix(NA, nrow = p, ncol = iter + 1L)
  g1_update[, 1] <- g1
  g2_update[, 1] <- g2

  AA <- integer(0)

  while (k < iter) {
    thetak <- thetaK

    # Compute dual norm and update direction
    step_info <- .compute_step(g1, g2, p, penalty, adap, eps2)
    delta <- step_info$delta
    i0 <- step_info$i0
    eps2_use <- step_info$eps2
    lambda_val <- step_info$lambda_val

    AA <- c(AA, i0)

    # Update theta
    thetaK <- thetak + (-1) * eps2_use * delta
    k <- k + 1L
    theta_update[, k + 1] <- thetaK

    alpha <- thetaK[1:p]         # recurrence
    beta  <- thetaK[(p + 1):(2 * p)]  # death

    # Update baseline hazards + recompute gradients in combined C++ calls
    res_de <- jfm_s0s1_cpp(Y, wt_matrix, td_sorted, index_death_matrix,
                            Z_pseudo, beta)
    S0t_de    <- res_de$S0t
    lambda0_d <- jfm_lambda0d_solution(td, d_td, n, S0t_de)

    res_re <- jfm_s0s1_cpp(Y, wt_recurrent_subject, tr_sorted,
                            index_recurrent_matrix, Z_pseudo, alpha)
    S0t_re    <- res_re$S0t
    lambda0_r <- jfm_lambda0r_solution(tr, d_tr, n, S0t_re)

    # Compute norm
    normK <- coop_norm(thetaK, p)
    normK_update[k + 1] <- normK

    # Scores from the same C++ call results
    g1 <- (-1) * jfm_score_cpp(index_recurrent_matrix, tr_id0,
                                 Z_pseudo, res_re$S1t, S0t_re) / n
    g2 <- (-1) * jfm_score_cpp(index_death_matrix, td_id0,
                                 Z_pseudo, res_de$S1t, S0t_de) / n

    # Gradient scaling: scale death (g2) up
    if (penalty %in% c("coop", "group")) {
      cc <- max(abs(g1)) / max(abs(g2))
    } else {
      cc <- sqrt(sum(g1^2)) / sqrt(sum(g2^2))
    }
    g2_origin <- g2
    g2 <- cc * g2

    norm2_g1_update[k + 1] <- sqrt(sum(g1^2))
    norm2_g2_update[k + 1] <- sqrt(sum(g2^2))
    g1_update[, k + 1] <- g1
    g2_update[, k + 1] <- g2

    lambda_vec[k] <- lambda_val

    # Early stop check: stop only if a single coordinate dominates every step
    # in the last pp iterations (truly stuck), not merely if two variables
    # alternate.
    if (k %% pp == 0) {
      if (length(unique(AA)) <= 1L) break
      AA <- integer(0)
    }
  }

  # Final lambda
  final_info <- .compute_step(g1, g2, p, penalty, adap, eps2)
  lambda_vec[k + 1] <- final_info$lambda_val

  # Trim storage
  idx_keep <- seq_len(k + 1)
  list(
    k = k,
    normK_update = normK_update[idx_keep],
    theta_update = theta_update[, idx_keep, drop = FALSE],
    norm2_g1_update = norm2_g1_update[idx_keep],
    lambda = lambda_vec[idx_keep],
    norm2_g2_update = norm2_g2_update[idx_keep],
    g1_update = g1_update[, idx_keep, drop = FALSE],
    g2_update = g2_update[, idx_keep, drop = FALSE]
  )
}


# --------------------------------------------------------------------------
# JSCM stagewise implementation
# --------------------------------------------------------------------------

#' @keywords internal
stagewise_jscm <- function(initial_alpha, initial_beta, Data2, penalty,
                           eps1, adap, eps2, iter, pp = 2000L) {
  p <- length(initial_alpha)
  n <- length(unique(Data2$id))

  # JSCM convention: alpha = recurrence (first p of theta),
  #                  beta = terminal/death (second p of theta)
  alpha <- initial_alpha
  beta <- initial_beta

  k <- 0L
  thetaK <- c(initial_alpha, initial_beta)
  normK <- 0

  # Initial gradients
  g1 <- (-1) * jscm_ee_alpha(Data2, alpha)
  g2 <- (-1) * jscm_ee_beta(Data2, alpha, beta)

  # Gradient scaling: scale death (g2) UP
  if (penalty %in% c("coop", "group")) {
    cc <- max(abs(g1)) / max(abs(g2))
  } else {
    cc <- sqrt(sum(g1^2)) / sqrt(sum(g2^2))
  }
  g2_origin <- g2
  g2 <- cc * g2

  # Storage
  theta_update <- matrix(NA, nrow = 2 * p, ncol = iter + 1L)
  theta_update[, 1] <- thetaK
  normK_update <- rep(NA_real_, iter + 1L)
  normK_update[1] <- normK
  norm2_g1_update <- rep(NA_real_, iter + 1L)
  norm2_g1_update[1] <- sqrt(sum(g1^2))
  norm2_g2_update <- rep(NA_real_, iter + 1L)
  norm2_g2_update[1] <- sqrt(sum(g2^2))
  lambda_vec <- rep(NA_real_, iter + 1L)
  g1_update <- matrix(NA, nrow = p, ncol = iter + 1L)
  g2_update <- matrix(NA, nrow = p, ncol = iter + 1L)
  g1_update[, 1] <- g1
  g2_update[, 1] <- g2

  AA <- integer(0)

  while (k < iter) {
    thetak <- thetaK

    # Compute dual norm and update direction
    step_info <- .compute_step(g1, g2, p, penalty, adap, eps2)
    delta <- step_info$delta
    i0 <- step_info$i0
    eps2_use <- step_info$eps2
    lambda_val <- step_info$lambda_val

    AA <- c(AA, i0)

    # Update theta
    thetaK <- thetak + (-1) * eps2_use * delta
    k <- k + 1L
    theta_update[, k + 1] <- thetaK

    alpha <- thetaK[1:p]
    beta <- thetaK[(p + 1):(2 * p)]

    # Compute norm
    normK <- coop_norm(thetaK, p)
    normK_update[k + 1] <- normK

    # Recompute gradients
    g1 <- (-1) * jscm_ee_alpha(Data2, alpha)
    g2 <- (-1) * jscm_ee_beta(Data2, alpha, beta)

    # Gradient scaling
    if (penalty %in% c("coop", "group")) {
      cc <- max(abs(g1)) / max(abs(g2))
    } else {
      cc <- sqrt(sum(g1^2)) / sqrt(sum(g2^2))
    }
    g2_origin <- g2
    g2 <- cc * g2

    norm2_g1_update[k + 1] <- sqrt(sum(g1^2))
    norm2_g2_update[k + 1] <- sqrt(sum(g2^2))
    g1_update[, k + 1] <- g1
    g2_update[, k + 1] <- g2

    lambda_vec[k] <- lambda_val

    # Early stop check: stop only if a single coordinate dominates every step
    # in the last pp iterations (truly stuck), not merely if two variables
    # alternate.
    if (k %% pp == 0) {
      if (length(unique(AA)) <= 1L) break
      AA <- integer(0)
    }
  }

  # Final lambda
  final_info <- .compute_step(g1, g2, p, penalty, adap, eps2)
  lambda_vec[k + 1] <- final_info$lambda_val

  # Trim storage
  idx_keep <- seq_len(k + 1)
  list(
    k = k,
    normK_update = normK_update[idx_keep],
    theta_update = theta_update[, idx_keep, drop = FALSE],
    norm2_g1_update = norm2_g1_update[idx_keep],
    lambda = lambda_vec[idx_keep],
    norm2_g2_update = norm2_g2_update[idx_keep],
    g1_update = g1_update[, idx_keep, drop = FALSE],
    g2_update = g2_update[, idx_keep, drop = FALSE]
  )
}


# --------------------------------------------------------------------------
# Shared helper: compute step direction and dual norm for all penalties
# --------------------------------------------------------------------------

#' @keywords internal
.compute_step <- function(g1, g2, p, penalty, adap, eps2) {
  delta <- numeric(2 * p)
  i0 <- 1L
  lambda_val <- 0
  eps2_use <- eps2

  if (penalty == "coop") {
    dual <- matrix(NA, nrow = p, ncol = 2)
    for (i in seq_len(p)) {
      if (g1[i] * g2[i] > 0) {
        dual[i, 1] <- sqrt(g1[i]^2 + g2[i]^2)
        dual[i, 2] <- 1
      } else {
        dual[i, 1] <- max(abs(g1[i]), abs(g2[i]))
        dual[i, 2] <- 2
      }
    }
    i0 <- which.max(dual[, 1])
    gi0 <- c(g1[i0], g2[i0])
    if (dual[i0, 2] == 1) {
      gi0 <- gi0 / sqrt(sum(gi0^2))
      delta[i0] <- gi0[1]
      delta[i0 + p] <- gi0[2]
    } else {
      if (abs(gi0[1]) >= abs(gi0[2])) {
        delta[i0] <- sign(gi0[1])
      } else {
        delta[p + i0] <- sign(gi0[2])
      }
    }
    lambda_val <- dual[i0, 1]
    if (adap == 1) {
      eps2_use <- 0.1^(count_digits(dual[i0, 1])) * 0.1
    }

  } else if (penalty == "lasso") {
    gg <- c(g1, g2)
    dual_val <- max(abs(gg))
    i0 <- which.max(abs(gg))
    delta[i0] <- sign(gg[i0])
    lambda_val <- dual_val
    if (adap == 1) {
      eps2_use <- 0.1^(count_digits(dual_val)) * 0.1
    }

  } else if (penalty == "group") {
    dual_vec <- numeric(p)
    for (i in seq_len(p)) {
      dual_vec[i] <- sqrt(g1[i]^2 + g2[i]^2)
    }
    i0 <- which.max(dual_vec)
    gi0 <- c(g1[i0], g2[i0])
    gi0 <- gi0 / sqrt(sum(gi0^2))
    delta[i0] <- gi0[1]
    delta[i0 + p] <- gi0[2]
    lambda_val <- dual_vec[i0]
    if (adap == 1) {
      eps2_use <- 0.1^(count_digits(dual_vec[i0])) * 0.1
    }
  }

  list(delta = delta, i0 = i0, eps2 = eps2_use, lambda_val = lambda_val)
}
