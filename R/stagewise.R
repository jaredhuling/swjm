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
#' @param eps Numeric. Base step size for the stagewise update, used at the top
#'   of the regularization path. The step is divided by ten each time the dual
#'   norm falls a further decade, so \code{eps} sets the resolution of the whole
#'   path. If \code{NULL}, defaults to 0.1 for the JFM and 0.01 for the JSCM.
#'   Smaller values trace the path more finely at proportionally greater cost;
#'   \code{max_iter} must be large enough for the path to reach small
#'   \code{lambda}.
#' @param max_iter Integer. Maximum number of stagewise iterations.
#' @param pp Integer. Early-stopping block size: algorithm checks every
#'   \code{pp} iterations if fewer than 3 unique coordinates were updated.
#' @param estimate_frailty Logical. For JFM only: if \code{TRUE}, estimates
#'   the frailty variance and uses the Kalbfleisch et al. (2013) frailty
#'   weights \eqn{w_i(t)} in the estimating equations. If \code{FALSE}
#'   (default), uses unit weights (simplified model without frailty).
#' @param lambda_min_ratio Numeric. The path stops once the dual norm falls
#'   below this fraction of its value at the top of the path. Past that point
#'   the coefficients are essentially static, so the remaining iterations cost
#'   time without changing the fit. If \code{NULL}, defaults to 0.01 when
#'   \code{n > p} and 1e-4 otherwise, following the convention used by
#'   \code{glmnet}. Set to 0 to trace the path to the smallest reachable
#'   \code{lambda}.
#' @param standardize Logical. If \code{TRUE} (default), covariates are
#'   divided by their standard deviations before fitting and coefficients
#'   are rescaled back to the original scale. For the JFM with time-varying
#'   covariates, the SD is computed across all rows (all time points and
#'   subjects). For the JSCM with time-invariant covariates, the SD is
#'   computed across subjects only.
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
                          eps = NULL, max_iter = NULL, pp = NULL,
                          estimate_frailty = FALSE,
                          standardize = TRUE,
                          lambda_min_ratio = NULL) {
  model <- match.arg(model)
  penalty <- match.arg(penalty)

  data <- prepare_data(data, caller = "stagewise_fit")
  validate_data(data, caller = "stagewise_fit")

  p <- ncol(data) - 5L
  cov_cols <- 6:ncol(data)
  initial_alpha <- numeric(p)
  initial_beta <- numeric(p)

  # Where to stop the path. Follows the glmnet convention: a short path when
  # there are many more observations than covariates, a longer one otherwise,
  # since with p >= n the interesting solutions sit further down the path.
  if (is.null(lambda_min_ratio)) {
    n_subj <- length(unique(data$id))
    lambda_min_ratio <- if (n_subj > p) 1e-2 else 1e-4
  }
  if (!is.numeric(lambda_min_ratio) || length(lambda_min_ratio) != 1L ||
      is.na(lambda_min_ratio) || lambda_min_ratio < 0 || lambda_min_ratio >= 1) {
    stop("'lambda_min_ratio' must be a single number in [0, 1).", call. = FALSE)
  }

  # Standardize covariates (divide by SD, no centering — Cox/AFT models
  # absorb the mean into the baseline)
  if (standardize) {
    if (model == "jfm") {
      # Time-varying: SD across ALL rows (all time points, all subjects)
      cov_sd <- apply(data[, cov_cols, drop = FALSE], 2, sd)
    } else {
      # Time-invariant (JSCM): SD across subjects (one terminal row per subject,
      # validated above)
      subj_rows <- !duplicated(data$id)
      cov_sd <- apply(data[subj_rows, cov_cols, drop = FALSE], 2, sd)
    }
    # Guard against zero-variance covariates
    cov_sd[cov_sd == 0] <- 1
    data_fit <- data
    data_fit[, cov_cols] <- sweep(data[, cov_cols, drop = FALSE], 2, cov_sd, "/")
  } else {
    data_fit <- data
    cov_sd <- rep(1, p)
  }

  if (model == "jfm") {
    if (is.null(eps)) eps <- 0.1
    if (is.null(max_iter)) max_iter <- 5000L
    if (is.null(pp)) pp <- max_iter   # disable early stopping by default for JFM
    result <- stagewise_jfm(initial_alpha, initial_beta, data_fit, penalty,
                            eps1 = 1e-6, adap = 1L, eps2 = eps,
                            iter = max_iter, pp = pp,
                            estimate_frailty = estimate_frailty,
                            lambda_min_ratio = lambda_min_ratio)
  } else {
    if (is.null(eps)) eps <- 0.01
    if (is.null(max_iter)) max_iter <- 5000L
    if (is.null(pp)) pp <- max_iter   # disable early stopping by default for JSCM
    result <- stagewise_jscm(initial_alpha, initial_beta, data_fit, penalty,
                             eps1 = 1e-6, adap = 1L, eps2 = eps,
                             iter = max_iter, pp = pp,
                             lambda_min_ratio = lambda_min_ratio)
  }

  # Rescale coefficients back to original covariate scale:
  # If z_std = z / sd, then alpha_std * z_std = (alpha_std / sd) * z,
  # so alpha_orig = alpha_std / sd.
  theta_raw <- result$theta_update
  theta_rescaled <- theta_raw
  theta_rescaled[1:p, ]         <- theta_raw[1:p, , drop = FALSE]         / cov_sd
  theta_rescaled[(p + 1):(2*p), ] <- theta_raw[(p + 1):(2*p), , drop = FALSE] / cov_sd

  structure(
    list(
      k = result$k,
      theta = theta_rescaled,
      alpha = theta_rescaled[1:p, , drop = FALSE],
      beta  = theta_rescaled[(p + 1):(2 * p), , drop = FALSE],
      lambda = result$lambda,
      norm = result$normK_update,
      norm_g1 = result$norm2_g1_update,
      norm_g2 = result$norm2_g2_update,
      g1 = result$g1_update,
      g2 = result$g2_update,
      model = model,
      penalty = penalty,
      p = p,
      cov_sd = cov_sd,
      standardize = standardize
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
                          eps1, adap, eps2, iter, pp = 200L,
                          estimate_frailty = FALSE, lambda_min_ratio = 0) {
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

  # Initialize baseline hazards and build pseudo-data via fast C++
  lambda0_d <- rep(0.2, length(td))
  theta <- 1
  Z1_raw <- do.call(rbind, Z)
  psd <- jfm_build_pseudo_cpp(t.start, I, Z1_raw, Y, td, td.id, tr, tr.id,
                                beta, lambda0_d, theta)

  Z_pseudo  <- psd$Z_pseudo
  td_sorted <- psd$td_sorted
  tr_sorted <- psd$tr_sorted
  td_id     <- as.integer(psd$td_id)
  tr_id     <- as.integer(psd$tr_id)
  de_epi    <- psd$de_epi
  re_epi    <- psd$re_epi
  entry_subject <- psd$entry_subject

  # Precompute timelines (one-time cost)
  exit_times  <- psd$exit_times
  is_last     <- psd$is_last
  entry_times <- psd$entry_times

  n_de <- length(td)
  n_re <- length(tr)

  if (estimate_frailty) {
    # Weighted timelines (include WEIGHT_UPDATE events at death times)
    tl_death <- jfm_precompute_timeline_wt(entry_times, exit_times, is_last,
                                            td_sorted, td_sorted, 0L)
    tl_recur <- jfm_precompute_timeline_wt(entry_times, exit_times, is_last,
                                            tr_sorted, td_sorted, 1L)
  } else {
    tl_death <- jfm_precompute_timeline(entry_times, exit_times, is_last,
                                         td_sorted, 0L)
    tl_recur <- jfm_precompute_timeline(entry_times, exit_times, is_last,
                                         tr_sorted, 1L)
  }

  # Helper: compute S0t/S1t + baseline hazards with current alpha, beta
  .compute_ee <- function(alpha, beta) {
    if (estimate_frailty) {
      # Update weights using current lambda0_d (warm start from previous step).
      # Two iterations suffice since coefficients change by only eps per step.
      for (wt_it in 1:2) {
        pp_wt <- jfm_build_pseudo_cpp(t.start, I, Z1_raw, Y,
          td, td.id, tr, tr.id, beta, lambda0_d, theta)
        wt_de <- pp_wt$wt_de
        r_de <- jfm_s0s1_wt_fast_cpp(tl_death$type, tl_death$idx, tl_death$size,
                                        Z_pseudo, entry_subject, beta, wt_de, n_de, n)
        lambda0_d <<- jfm_lambda0d_solution(td, d_td, n, r_de$S0t)
      }
      r_de$score <- jfm_score_fast_cpp(de_epi, Z_pseudo, r_de$S1t, r_de$S0t) / n
      r_re <- jfm_s0s1_wt_fast_cpp(tl_recur$type, tl_recur$idx, tl_recur$size,
                                      Z_pseudo, entry_subject, alpha, wt_de, n_re, n)
      r_re$score <- jfm_score_fast_cpp(re_epi, Z_pseudo, r_re$S1t, r_re$S0t) / n
    } else {
      r_de <- jfm_s0s1_score_fused_cpp(tl_death$type, tl_death$idx, tl_death$size,
                                          Z_pseudo, beta, de_epi, n_de, n)
      r_re <- jfm_s0s1_score_fused_cpp(tl_recur$type, tl_recur$idx, tl_recur$size,
                                          Z_pseudo, alpha, re_epi, n_re, n)
    }
    list(res_de = r_de, res_re = r_re)
  }

  # Initial computation
  ee <- .compute_ee(alpha, beta)
  S0t_de <- ee$res_de$S0t
  lambda0_d <- jfm_lambda0d_solution(td, d_td, n, S0t_de)
  S0t_re <- ee$res_re$S0t
  lambda0_r <- jfm_lambda0r_solution(tr, d_tr, n, S0t_re)

  # Stagewise loop
  k <- 0L
  thetaK <- c(initial_alpha, initial_beta)
  normK <- 0

  # Initial gradients (fused score from .compute_ee)
  g1 <- (-1) * as.numeric(ee$res_re$score)
  g2 <- (-1) * as.numeric(ee$res_de$score)

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
  lam0 <- NA_real_

  while (k < iter) {
    thetak <- thetaK

    # Compute dual norm and update direction. The first pass fixes lam0, the
    # dual norm at the top of the path, which anchors the adaptive step size.
    step_info <- .compute_step(g1, g2, p, penalty, adap, eps2, cc, lam0)
    if (!is.finite(lam0)) {
      lam0 <- step_info$lambda_val
      step_info <- .compute_step(g1, g2, p, penalty, adap, eps2, cc, lam0)
    }
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

    # Update baseline hazards + recompute gradients
    ee <- .compute_ee(alpha, beta)
    S0t_de    <- ee$res_de$S0t
    lambda0_d <- jfm_lambda0d_solution(td, d_td, n, S0t_de)
    S0t_re    <- ee$res_re$S0t
    lambda0_r <- jfm_lambda0r_solution(tr, d_tr, n, S0t_re)

    # Compute norm
    normK <- coop_norm(thetaK, p)
    normK_update[k + 1] <- normK

    # Scores (fused with S0t computation)
    g1 <- (-1) * as.numeric(ee$res_re$score)
    g2 <- (-1) * as.numeric(ee$res_de$score)

    # Rescale the death gradient with cc FROZEN at its initial value, so the
    # algorithm emulates a single fixed scaled penalty (recalibrating cc every
    # iteration would feed overshoot back into the amplified steps and
    # destabilize the path)
    g2_origin <- g2
    g2 <- cc * g2

    norm2_g1_update[k + 1] <- sqrt(sum(g1^2))
    norm2_g2_update[k + 1] <- sqrt(sum(g2^2))
    g1_update[, k + 1] <- g1
    g2_update[, k + 1] <- g2

    lambda_vec[k] <- lambda_val

    # Stop once the dual norm has fallen to lambda_min_ratio of its value at
    # the top of the path. Beyond that point the coefficients are essentially
    # static -- the remaining iterations contribute a negligible fraction of
    # the total movement -- so continuing only costs time and inflates the
    # lambda grid that cross-validation must search.
    if (is.finite(lam0) && lambda_val < lambda_min_ratio * lam0) break

    # Early stop check: stop only if a single coordinate dominates every step
    # in the last pp iterations (truly stuck), not merely if two variables
    # alternate.
    if (k %% pp == 0) {
      if (length(unique(AA)) <= 1L) break
      AA <- integer(0)
    }
  }

  # Final lambda
  final_info <- .compute_step(g1, g2, p, penalty, adap, eps2, cc, lam0)
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
                           eps1, adap, eps2, iter, pp = 2000L,
                           lambda_min_ratio = 0) {
  p <- length(initial_alpha)
  n <- length(unique(Data2$id))

  # JSCM convention: alpha = recurrence (first p of theta),
  #                  beta = terminal/death (second p of theta)
  alpha <- initial_alpha
  beta <- initial_beta

  # --- Precompute data structures once (avoid recomputing each iteration) ---
  xi <- as.matrix(Data2[Data2$event == 0, 6:ncol(Data2)])
  rownames(xi) <- seq_len(nrow(xi))
  ti <- Data2$t.stop[Data2$event == 1]
  yi <- Data2$t.stop[Data2$event == 0]
  di <- Data2$status[Data2$event == 0]

  uids <- unique(Data2$id)
  Data_recur <- Data2[Data2$event == 1, ]
  xi_recur <- as.matrix(Data_recur[, 6:ncol(Data2)])
  m <- as.integer(vapply(uids, function(id) sum(Data_recur$id == id), integer(1)))
  w1 <- rep(1, length(yi))

  .jscm_g1 <- function(a) {
    out <- am1(a, ti, yi, w1, xi, m)
    (-1) * (if (is.matrix(out)) drop(out) else as.numeric(out))
  }

  .jscm_g2 <- function(a, b) {
    texa <- log(ti) + xi_recur %*% a
    yexa <- log(yi) + xi %*% a
    rate <- c(reRate(texa, rep(yexa, m), rep(w1, m), yexa))
    Lam <- exp(-rate)
    numAdj <- 1e-3
    if (numAdj > min(Lam)) numAdj <- numAdj * min(Lam)
    R2 <- (m + numAdj) / (Lam + numAdj)
    out <- temLog(b, b, xi, yi, R2, di, w1)
    (-1) * (if (is.matrix(out)) drop(out) else as.numeric(out))
  }

  k <- 0L
  thetaK <- c(initial_alpha, initial_beta)
  normK <- 0

  # Initial gradients
  g1 <- .jscm_g1(alpha)
  g2 <- .jscm_g2(alpha, beta)

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
  lam0 <- NA_real_

  while (k < iter) {
    thetak <- thetaK

    # Compute dual norm and update direction. The first pass fixes lam0, the
    # dual norm at the top of the path, which anchors the adaptive step size.
    step_info <- .compute_step(g1, g2, p, penalty, adap, eps2, cc, lam0)
    if (!is.finite(lam0)) {
      lam0 <- step_info$lambda_val
      step_info <- .compute_step(g1, g2, p, penalty, adap, eps2, cc, lam0)
    }
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

    # Recompute gradients (using precomputed data structures)
    g1 <- .jscm_g1(alpha)
    g2 <- .jscm_g2(alpha, beta)

    # Rescale the death gradient with cc FROZEN at its initial value, so the
    # algorithm emulates a single fixed scaled penalty (recalibrating cc every
    # iteration would feed overshoot back into the amplified steps and
    # destabilize the path)
    g2_origin <- g2
    g2 <- cc * g2

    norm2_g1_update[k + 1] <- sqrt(sum(g1^2))
    norm2_g2_update[k + 1] <- sqrt(sum(g2^2))
    g1_update[, k + 1] <- g1
    g2_update[, k + 1] <- g2

    lambda_vec[k] <- lambda_val

    # Stop once the dual norm has fallen to lambda_min_ratio of its value at
    # the top of the path. Beyond that point the coefficients are essentially
    # static -- the remaining iterations contribute a negligible fraction of
    # the total movement -- so continuing only costs time and inflates the
    # lambda grid that cross-validation must search.
    if (is.finite(lam0) && lambda_val < lambda_min_ratio * lam0) break

    # Early stop check: stop only if a single coordinate dominates every step
    # in the last pp iterations (truly stuck), not merely if two variables
    # alternate.
    if (k %% pp == 0) {
      if (length(unique(AA)) <= 1L) break
      AA <- integer(0)
    }
  }

  # Final lambda
  final_info <- .compute_step(g1, g2, p, penalty, adap, eps2, cc, lam0)
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
#'
#' g1 and g2 are the (loss-convention) gradients with g2 already rescaled by
#' cc = 1/xi, so coordinate selection and the dual-norm value are computed on
#' the rescaled pair. Under the scaled penalties, the exact argmax direction
#' carries an additional factor 1/xi = cc on every death (second-block)
#' coordinate: the scaled lasso and the mixed-sign cooperative branch give
#' e_{j*} sgn(g)/xi when a death coordinate wins, and the group / same-sign
#' cooperative L2 direction is (g_a, g_b/xi^2)/||(g_a, g_b/xi)||_2, i.e. cc
#' times the rescaled component. The amplification is applied at the end.
#' @keywords internal
#'
#' Adaptive step size. The step starts at \code{eps2} at the top of the path
#' and is divided by ten each time the dual norm falls a further decade below
#' its initial value \code{lam0}. Anchoring to \code{lam0} (rather than to the
#' absolute magnitude of the dual norm) keeps the rule scale invariant, so the
#' step depends only on how far along the path we are, and makes \code{eps}
#' the meaningful base resolution for both model types.
.adaptive_eps <- function(eps2, dual, lam0) {
  if (!is.finite(lam0) || lam0 <= 0 || !is.finite(dual) || dual <= 0) return(eps2)
  k <- max(0, floor(log10(lam0)) - floor(log10(dual)))
  eps2 * 10^(-k)
}

.compute_step <- function(g1, g2, p, penalty, adap, eps2, cc = 1, lam0 = NA_real_) {
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
      eps2_use <- .adaptive_eps(eps2, dual[i0, 1], lam0)
    }

  } else if (penalty == "lasso") {
    gg <- c(g1, g2)
    dual_val <- max(abs(gg))
    i0 <- which.max(abs(gg))
    delta[i0] <- sign(gg[i0])
    lambda_val <- dual_val
    if (adap == 1) {
      eps2_use <- .adaptive_eps(eps2, dual_val, lam0)
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
      eps2_use <- .adaptive_eps(eps2, dual_vec[i0], lam0)
    }
  }

  if (cc != 1) {
    delta[(p + 1):(2 * p)] <- cc * delta[(p + 1):(2 * p)]
  }

  list(delta = delta, i0 = i0, eps2 = eps2_use, lambda_val = lambda_val)
}
