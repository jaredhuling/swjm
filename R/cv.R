# ============================================================================
# Cross-validation for both JFM and JSCM
# ============================================================================


#' Cross-Validation for Stagewise Variable Selection
#'
#' Selects the optimal penalty parameter (lambda) along the stagewise path
#' using K-fold cross-validation with cross-fitted estimating equations.
#'
#' @param data A data frame in recurrent-event format.
#' @param model Character. Either \code{"jfm"} or \code{"jscm"}.
#' @param penalty Character. One of \code{"coop"}, \code{"lasso"}, or
#'   \code{"group"}.
#' @param K Integer. Number of cross-validation folds (default 5).
#' @param lambda_seq Numeric vector. The lambda sequence at which to
#'   evaluate the cross-validation criterion. If \code{NULL}, extracted
#'   from a full-data stagewise fit.
#' @param eps Numeric. Step size (passed to \code{stagewise_fit}).
#' @param max_iter Integer. Maximum iterations (passed to \code{stagewise_fit}).
#' @param pp Integer. Early-stop block size (passed to \code{stagewise_fit}).
#'
#' @return An object of class \code{"swjm_cv"}, a list with components:
#'   \describe{
#'     \item{position_CF}{Integer, position of best lambda by combined
#'       cross-fitted score norm.}
#'     \item{position_CF_re}{Integer, position of best lambda by
#'       recurrence score norm.}
#'     \item{position_CF_cen}{Integer, position of best lambda by
#'       terminal score norm.}
#'     \item{lambda_seq}{Numeric vector of lambda values evaluated.}
#'     \item{Scorenorm_crossfit}{Combined cross-fitted score norm path.}
#'     \item{Scorenorm_crossfit_re}{Recurrence score norm path.}
#'     \item{Scorenorm_crossfit_ce}{Terminal score norm path.}
#'     \item{full_fit}{The full-data stagewise fit (class \code{"swjm_path"}).}
#'     \item{model}{Character.}
#'     \item{penalty}{Character.}
#'   }
#'
#' @examples
#' \donttest{
#' dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
#' fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
#'                      max_iter = 100)
#' cv_res <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
#'                        lambda_seq = fit$lambda, max_iter = 100)
#' cv_res
#' }
#'
#' @export
cv_stagewise <- function(data, model = c("jfm", "jscm"),
                         penalty = c("coop", "lasso", "group"),
                         K = 5L, lambda_seq = NULL,
                         eps = NULL, max_iter = NULL, pp = NULL) {
  model <- match.arg(model)
  penalty <- match.arg(penalty)

  p <- ncol(data) - 5L
  initial_alpha <- numeric(p)
  initial_beta <- numeric(p)

  # Set defaults
  if (model == "jfm") {
    if (is.null(eps)) eps <- 0.1
    if (is.null(max_iter)) max_iter <- 5000L
    if (is.null(pp)) pp <- max_iter   # disable early stopping by default for JFM
  } else {
    if (is.null(eps)) eps <- 0.01
    if (is.null(max_iter)) max_iter <- 5000L
    if (is.null(pp)) pp <- max_iter   # disable early stopping by default for JSCM
  }

  # Always run full-data fit (needed for coef extraction)
  full_fit <- stagewise_fit(data, model = model, penalty = penalty,
                            eps = eps, max_iter = max_iter, pp = pp)

  if (is.null(lambda_seq)) {
    lambda_seq <- full_fit$lambda
    # Extract decreasing path
    dec_idx <- extract_decreasing_indices(lambda_seq)
    lambda_seq <- lambda_seq[dec_idx]
  }

  if (model == "jfm") {
    result <- cv_jfm(data, penalty, lambda_seq, K, initial_alpha,
                     initial_beta, eps1 = 1e-6, adap = 1L, eps2 = eps,
                     iter = max_iter, pp = pp)
  } else {
    result <- cv_jscm(data, penalty, lambda_seq, K, initial_alpha,
                      initial_beta, eps1 = 1e-6, adap = 1L, eps2 = eps,
                      iter = max_iter, pp = pp)
  }

  # --- Derived quantities from the full-data path ---
  lambda_path <- full_fit$lambda
  dec_idx     <- extract_decreasing_indices(lambda_path)
  theta_dec   <- full_fit$theta[, dec_idx, drop = FALSE]
  lambda_dec  <- lambda_path[dec_idx]

  # Best alpha / beta at optimal lambda
  target_lambda <- lambda_seq[result$position_CF]
  lam_best      <- lambda_interp(lambda_dec, target_lambda)
  theta_best    <- as.numeric(theta_dec[, lam_best$left])  * lam_best$frac +
                   as.numeric(theta_dec[, lam_best$right]) * (1 - lam_best$frac)
  alpha_best <- theta_best[1:p]
  beta_best  <- theta_best[(p + 1):(2 * p)]

  # Variable counts along the full lambda_seq
  lam_all   <- lambda_interp(lambda_dec, lambda_seq)
  theta_all <- as.matrix(
    theta_dec[, lam_all$left,  drop = FALSE] %*% Matrix::Diagonal(x = lam_all$frac) +
    theta_dec[, lam_all$right, drop = FALSE] %*% Matrix::Diagonal(x = 1 - lam_all$frac)
  )
  n_active_alpha <- colSums(theta_all[1:p, , drop = FALSE]             != 0)
  n_active_beta  <- colSums(theta_all[(p + 1):(2 * p), , drop = FALSE] != 0)
  n_active       <- n_active_alpha + n_active_beta

  # Baseline hazards — JFM via Breslow, JSCM via Nelson-Aalen on accelerated scale
  if (model == "jfm") {
    dc_base  <- extract_data_components(data)
    haz      <- lambda0_gen(data, alpha_best, beta_best)
    baseline <- list(
      lambda0_r = haz$lambda0_r,
      lambda0_d = haz$lambda0_d,
      tr = sort(dc_base$tr),
      td = sort(dc_base$td)
    )
  } else {
    baseline <- .jscm_baseline(data, alpha_best, beta_best)
  }

  structure(
    c(result,
      list(
        lambda_seq     = lambda_seq,
        full_fit       = full_fit,
        model          = model,
        penalty        = penalty,
        p              = p,
        alpha          = alpha_best,
        beta           = beta_best,
        n_active_alpha = n_active_alpha,
        n_active_beta  = n_active_beta,
        n_active       = n_active,
        baseline       = baseline
      )),
    class = "swjm_cv"
  )
}


#' @export
print.swjm_cv <- function(x, ...) {
  cat("Cross-validation (", x$model, "/", x$penalty, ")\n\n", sep = "")
  cat(sprintf("  %-28s %d\n", "Covariates (p):", x$p))
  cat(sprintf("  %-28s %d\n", "Lambda grid size:", length(x$lambda_seq)))
  cat(sprintf("  %-28s %d  (lambda = %.4g)\n", "Best position (combined):",
              x$position_CF, x$lambda_seq[x$position_CF]))
  na <- sum(x$alpha != 0)
  nb <- sum(x$beta  != 0)
  cat(sprintf("  %-28s %d readmission, %d death\n",
              "Selected variables:", na, nb))
  if (na > 0)
    cat(sprintf("    Readmission (alpha): %s\n",
                paste(which(x$alpha != 0), collapse = ", ")))
  if (nb > 0)
    cat(sprintf("    Death (beta):        %s\n",
                paste(which(x$beta != 0), collapse = ", ")))
  invisible(x)
}


#' @export
coef.swjm_cv <- function(object, ...) {
  c(object$alpha, object$beta)
}


#' @export
summary.swjm_cv <- function(object, ...) {
  p <- object$p

  cat(sprintf("CV-selected model (%s/%s)\n\n", object$model, object$penalty))
  cat(sprintf("  p = %d  |  Lambda grid: %d steps  |  CV optimal: step %d (lambda = %.4g)\n\n",
              p, length(object$lambda_seq),
              object$position_CF, object$lambda_seq[object$position_CF]))

  a <- object$alpha
  b <- object$beta
  sel <- which(a != 0 | b != 0)

  if (length(sel) == 0L) {
    cat("  No variables selected.\n")
    return(invisible(object))
  }

  # Sort by decreasing |alpha| + |beta|
  ord <- sel[order(abs(a[sel]) + abs(b[sel]), decreasing = TRUE)]

  .var_type <- function(j) {
    if (a[j] != 0 && b[j] != 0) {
      if (sign(a[j]) == sign(b[j])) "shared (+)" else "shared (\u2013)"
    } else if (a[j] != 0) {
      "readmission only"
    } else {
      "death only"
    }
  }

  cat(sprintf("  Selected coefficients  (%d readmission, %d death):\n\n",
              sum(a != 0), sum(b != 0)))
  cat(sprintf("  %-10s  %-10s  %-10s  %s\n",
              "Variable", "alpha", "beta", "Type"))
  cat("  ----------  ----------  ----------  ----------------\n")
  for (j in ord) {
    a_str <- if (a[j] != 0) sprintf("%+.4f", a[j]) else "     \u2014"
    b_str <- if (b[j] != 0) sprintf("%+.4f", b[j]) else "     \u2014"
    cat(sprintf("  %-10s  %-10s  %-10s  %s\n",
                paste0("x", j), a_str, b_str, .var_type(j)))
  }

  inactive <- setdiff(seq_len(p), sel)
  if (length(inactive) > 0L)
    cat(sprintf("\n  Inactive (%d): %s\n", length(inactive),
                paste(paste0("x", inactive), collapse = ", ")))

  invisible(object)
}


# --------------------------------------------------------------------------
# JFM cross-validation
# --------------------------------------------------------------------------

#' @keywords internal
cv_jfm <- function(Data2, penalty, lambda_seq, K, initial_alpha,
                   initial_beta, eps1, adap, eps2, iter, pp) {
  p <- ncol(Data2) - 5L
  n1 <- length(unique(Data2$id))

  # Stratified CV folds
  folds <- .stratified_folds(Data2, K)

  # Store fold-wise solution paths and fold indices for cross-fitting
  test_inx <- matrix(NA_integer_, nrow = n1, ncol = 2)
  theta_lambda_list <- vector("list", K)

  for (kk in seq_len(K)) {
    leave_out <- folds[[kk]]$test_ids
    test_inx[leave_out, 1] <- kk
    test_inx[leave_out, 2] <- leave_out

    Data2_tr <- Data2[!(Data2$id %in% leave_out), ]
    rownames(Data2_tr) <- seq_len(nrow(Data2_tr))
    Data2_tr$id <- match(Data2_tr$id, unique(Data2_tr$id))

    # Train on fold
    results_tr <- stagewise_jfm(initial_alpha, initial_beta, Data2_tr,
                                penalty, eps1, adap, eps2, iter, pp)

    lambda_seq_tr <- results_tr$lambda
    theta_tr <- results_tr$theta_update

    # Extract decreasing lambda path
    dec_idx <- extract_decreasing_indices(lambda_seq_tr)
    lambda_seq_tr_ordered <- lambda_seq_tr[dec_idx]
    theta_tr_ordered <- theta_tr[, dec_idx, drop = FALSE]

    # Interpolate to evaluation lambda_seq
    lamlist <- lambda_interp(lambda_seq_tr_ordered, lambda_seq)
    theta_lambda <- as.matrix(
      theta_tr_ordered[, lamlist$left, drop = FALSE] %*%
        Matrix::Diagonal(x = lamlist$frac) +
        theta_tr_ordered[, lamlist$right, drop = FALSE] %*%
        Matrix::Diagonal(x = 1 - lamlist$frac)
    )
    theta_lambda_list[[kk]] <- theta_lambda
  }

  # Cross-fitted EE evaluation on full data
  dc <- extract_data_components(Data2)
  Z <- dc$Z; n <- dc$n; td <- dc$td; td.id <- dc$td.id; d_td <- dc$d_td
  tr <- dc$tr; tr.id <- dc$tr.id; d_tr <- dc$d_tr
  Y <- dc$Y; STATUS <- dc$STATUS
  t.start <- dc$t.start; I <- dc$I

  lambda0_d_init <- rep(0.2, length(td))
  theta_param <- 1
  alpha0 <- beta0 <- numeric(p)

  A <- jfm_wt_death(theta_param, alpha0, t.start, I, Z, td,
                     lambda0_d_init, td.id)
  td_id              <- A$td_id
  index_death_matrix <- A$index_death_matrix
  pseudo_entries     <- A$pseudo_entries
  wt_matrix          <- matrix(1, nrow = n, ncol = length(td))

  diff_tr1       <- diff(c(0, tr[order(tr)]))
  lambda0_r_init <- diff_tr1 * 1
  B <- jfm_r2i_integral(t.start, I, Z, beta0, tr, lambda0_r_init, tr.id)
  index_recurrent_matrix <- B$index_recurrent_matrix
  tr_id                  <- B$tr_id
  wt_recurrent_subject   <- matrix(1, nrow = n, ncol = length(tr))

  # Pre-extract shared objects used in the hot CV loop
  Z_pseudo  <- pseudo_entries[, 3:ncol(pseudo_entries), drop = FALSE]
  td_sorted <- sort(td)
  tr_sorted <- sort(tr)
  td_id0    <- td_id - 1L   # 0-based for C++ score
  tr_id0    <- tr_id - 1L
  # cv_fold: 0-based fold assignment indexed by subject ID (test_inx[subj, 1])
  cv_fold   <- as.integer(test_inx[, 1]) - 1L   # 0-based

  # Scaling factor from null model (C++ combined call)
  res_de0 <- jfm_s0s1_cpp(Y, wt_matrix, td_sorted, index_death_matrix,
                            Z_pseudo, beta0)
  res_re0 <- jfm_s0s1_cpp(Y, wt_recurrent_subject, tr_sorted,
                            index_recurrent_matrix, Z_pseudo, alpha0)
  g10 <- (-1) * jfm_score_cpp(index_recurrent_matrix, tr_id0,
                                Z_pseudo, res_re0$S1t, res_re0$S0t) / n
  g20 <- (-1) * jfm_score_cpp(index_death_matrix, td_id0,
                                Z_pseudo, res_de0$S1t, res_de0$S0t) / n
  c1 <- max(abs(g10)) / max(abs(g20))

  # Cross-fitted score norms
  Scorenorm_crossfit_re  <- numeric(length(lambda_seq))
  Scorenorm_crossfit_cen <- numeric(length(lambda_seq))
  Scorenorm_crossfit     <- numeric(length(lambda_seq))

  for (j in seq_along(lambda_seq)) {
    alpha_mat <- matrix(NA, nrow = K, ncol = p)
    beta_mat  <- matrix(NA, nrow = K, ncol = p)
    for (qq in seq_len(K)) {
      alpha_mat[qq, ] <- theta_lambda_list[[qq]][1:p, j]
      beta_mat[qq, ]  <- theta_lambda_list[[qq]][(p + 1):(2 * p), j]
    }

    # C++ cross-fitted S0t+S1t (one call replaces four R loops each)
    res_de_cf <- jfm_s0s1_cf_cpp(Y, wt_matrix, td_sorted,
                                   index_death_matrix, Z_pseudo,
                                   beta_mat, cv_fold)
    res_re_cf <- jfm_s0s1_cf_cpp(Y, wt_recurrent_subject, tr_sorted,
                                   index_recurrent_matrix, Z_pseudo,
                                   alpha_mat, cv_fold)

    g1 <- (-1) * jfm_score_cpp(index_recurrent_matrix, tr_id0,
                                 Z_pseudo, res_re_cf$S1t, res_re_cf$S0t) / n
    g2 <- (-1) * jfm_score_cpp(index_death_matrix, td_id0,
                                 Z_pseudo, res_de_cf$S1t, res_de_cf$S0t) / n
    g2_origin <- g2
    g2 <- g2 * c1

    Scorenorm_crossfit_re[j]  <- sqrt(sum(g1^2))
    Scorenorm_crossfit_cen[j] <- sqrt(sum(g2_origin^2))
    Scorenorm_crossfit[j]     <- sqrt(sum(c(g1, g2)^2))
  }

  list(
    position_CF_re = which.min(Scorenorm_crossfit_re),
    position_CF_cen = which.min(Scorenorm_crossfit_cen),
    position_CF = which.min(Scorenorm_crossfit),
    Scorenorm_crossfit_re = Scorenorm_crossfit_re,
    Scorenorm_crossfit_ce = Scorenorm_crossfit_cen,
    Scorenorm_crossfit = Scorenorm_crossfit
  )
}


# --------------------------------------------------------------------------
# JSCM cross-validation
# --------------------------------------------------------------------------

#' @keywords internal
cv_jscm <- function(Data2, penalty, lambda_seq, K, initial_alpha,
                    initial_beta, eps1, adap, eps2, iter, pp) {
  p <- ncol(Data2) - 5L
  n1 <- length(unique(Data2$id))

  # Stratified CV folds
  folds <- .stratified_folds(Data2, K)

  # Collect info from each fold for cross-fitting EE
  num_re <- sum(Data2$event)
  ti_vec <- numeric(num_re)
  w1_vec <- numeric(n1)
  xi_mat <- matrix(NA, nrow = n1, ncol = p)
  m_vec <- numeric(n1)
  di_vec <- numeric(n1)

  texa_mat <- matrix(NA, nrow = num_re, ncol = length(lambda_seq))
  yexa_mat <- matrix(NA, nrow = n1, ncol = length(lambda_seq))
  yexb_mat <- matrix(NA, nrow = n1, ncol = length(lambda_seq))

  count_re <- 0L
  count_cen <- 0L

  for (kk in seq_len(K)) {
    leave_out <- folds[[kk]]$test_ids

    Data2_tr <- Data2[!(Data2$id %in% leave_out), ]
    rownames(Data2_tr) <- seq_len(nrow(Data2_tr))
    Data2_tr$id <- match(Data2_tr$id, unique(Data2_tr$id))

    Data2_val <- Data2[Data2$id %in% leave_out, ]
    rownames(Data2_val) <- seq_len(nrow(Data2_val))
    Data2_val$id <- match(Data2_val$id, unique(Data2_val$id))

    # Train
    results_tr <- stagewise_jscm(initial_alpha, initial_beta, Data2_tr,
                                 penalty, eps1, adap, eps2, iter, pp)
    lambda_seq_tr <- results_tr$lambda
    theta_tr <- results_tr$theta_update

    # Extract decreasing path
    dec_idx <- extract_decreasing_indices(lambda_seq_tr)
    lambda_seq_tr_ordered <- lambda_seq_tr[dec_idx]
    theta_tr_ordered <- theta_tr[, dec_idx, drop = FALSE]

    # Interpolate
    lamlist <- lambda_interp(lambda_seq_tr_ordered, lambda_seq)
    theta_lambda <- as.matrix(
      theta_tr_ordered[, lamlist$left, drop = FALSE] %*%
        Matrix::Diagonal(x = lamlist$frac) +
        theta_tr_ordered[, lamlist$right, drop = FALSE] %*%
        Matrix::Diagonal(x = 1 - lamlist$frac)
    )

    # Collect test set info
    count_re_pre <- count_re
    count_re <- count_re + sum(Data2_val$event)
    count_cen_pre <- count_cen
    n_val <- sum(Data2_val$event == 0)
    count_cen <- count_cen + n_val

    ti <- Data2_val$t.stop[Data2_val$event == 1]
    ti_vec[(count_re_pre + 1):count_re] <- ti
    w1_vec[(count_cen_pre + 1):count_cen] <- rep(1, n_val)
    xi <- as.matrix(Data2_val[Data2_val$event == 0, 6:ncol(Data2_val)])
    xi_mat[(count_cen_pre + 1):count_cen, ] <- xi

    uids_val <- unique(Data2_val$id)
    Data_recur <- Data2_val[Data2_val$event == 1, ]
    list_recur_val <- vector("list", n_val)
    for (i2 in seq_len(n_val)) {
      list_recur_val[[i2]] <- Data_recur$t.stop[Data_recur$id == uids_val[i2]]
    }
    m <- as.integer(vapply(list_recur_val, length, integer(1)))
    m_vec[(count_cen_pre + 1):count_cen] <- m
    di_vec[(count_cen_pre + 1):count_cen] <- Data2_val$status[Data2_val$event == 0]

    for (j in seq_along(lambda_seq)) {
      alpha_j <- theta_lambda[1:p, j]
      beta_j <- theta_lambda[(p + 1):(2 * p), j]
      yi <- Data2_val$t.stop[Data2_val$event == 0]
      texa <- drop(log(ti) + as.matrix(Data2_val[Data2_val$event == 1, -(1:5)]) %*% alpha_j)
      yexa <- drop(log(yi) + xi %*% alpha_j)
      yexb <- drop(yi * exp(xi %*% beta_j))

      texa_mat[(count_re_pre + 1):count_re, j] <- texa
      yexa_mat[(count_cen_pre + 1):count_cen, j] <- yexa
      yexb_mat[(count_cen_pre + 1):count_cen, j] <- yexb
    }
  }

  # Cross-fitted EE using pooled test-set info
  Scorenorm_crossfit_re <- numeric(length(lambda_seq))
  Scorenorm_crossfit_cen <- numeric(length(lambda_seq))
  Scorenorm_crossfit <- numeric(length(lambda_seq))

  Tvec <- ti_vec
  W <- w1_vec
  X <- xi_mat
  m_all <- m_vec

  # Scaling factor from null model
  g10 <- (-1) * jscm_ee_alpha(Data2, numeric(p))
  g20 <- (-1) * jscm_ee_beta(Data2, numeric(p), numeric(p))
  c1 <- max(abs(g10)) / max(abs(g20))

  for (j in seq_along(lambda_seq)) {
    # EE of alpha (recurrence)
    texa <- texa_mat[, j]
    yexa <- yexa_mat[, j]
    out <- am2(texa, rep(yexa, m_all), yexa, Tvec, W, X, m_all)
    g1 <- (-1) * (if (is.matrix(out)) drop(out) else as.numeric(out))

    # EE of beta (death)
    rate <- c(reRate(texa, rep(yexa, m_all), rep(W, m_all), yexa))
    Lam <- exp(-rate)
    R <- m_all / Lam
    numAdj <- 1e-3
    if (numAdj > min(Lam)) numAdj <- numAdj * min(Lam)
    R2 <- (m_all + numAdj) / (Lam + numAdj)

    yexb <- yexb_mat[, j]
    di <- di_vec
    out2 <- temLog2(yexb, X, R2, di, W)
    g2 <- if (is.matrix(out2)) drop(out2) else as.numeric(out2)
    g2_origin <- g2
    g2 <- g2 * c1 * (-1)

    Scorenorm_crossfit_re[j] <- sqrt(sum(g1^2))
    Scorenorm_crossfit_cen[j] <- sqrt(sum(g2_origin^2))
    Scorenorm_crossfit[j] <- sqrt(sum(c(g1, g2)^2))
  }

  list(
    position_CF_re = which.min(Scorenorm_crossfit_re),
    position_CF_cen = which.min(Scorenorm_crossfit_cen),
    position_CF = which.min(Scorenorm_crossfit),
    Scorenorm_crossfit_re = Scorenorm_crossfit_re,
    Scorenorm_crossfit_ce = Scorenorm_crossfit_cen,
    Scorenorm_crossfit = Scorenorm_crossfit
  )
}


# --------------------------------------------------------------------------
# Shared helper: stratified CV folds
# --------------------------------------------------------------------------

#' @keywords internal
.stratified_folds <- function(Data2, K) {
  id_status <- stats::aggregate(status ~ id, data = Data2, FUN = sum)
  id_status$class <- as.integer(id_status$status > 0)

  id_class0 <- id_status$id[id_status$class == 0]
  id_class1 <- id_status$id[id_status$class == 1]

  folds0 <- create_folds(id_class0, K)
  folds1 <- create_folds(id_class1, K)

  lapply(seq_len(K), function(kk) {
    test_ids <- c(folds0[[kk]], folds1[[kk]])
    list(test_ids = test_ids)
  })
}
