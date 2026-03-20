# ============================================================================
# predict.swjm_cv and associated plot method
# ============================================================================


#' Predict from a Cross-Validated Joint Model Fit
#'
#' Computes subject-specific predictions from a cross-validated \code{swjm_cv}
#' fit.  The output differs by model:
#'
#' \strong{JFM} — returns survival curves for both processes using the Breslow
#' cumulative baseline hazards.  For subject \eqn{i} with covariate vector
#' \eqn{z_i}:
#' \deqn{
#'   S_{\text{re}}(t \mid z_i) = \exp\!\bigl(-\hat\Lambda_0^r(t)\,
#'     e^{\hat\alpha^\top z_i}\bigr), \quad
#'   S_{\text{de}}(t \mid z_i) = \exp\!\bigl(-\hat\Lambda_0^d(t)\,
#'     e^{\hat\beta^\top z_i}\bigr).
#' }
#'
#' \strong{JSCM} — returns survival curves for both processes using a
#' Nelson-Aalen baseline on the accelerated time scale.  For subject \eqn{i}:
#' \deqn{
#'   S_{\text{re}}(t \mid z_i) = \exp\!\bigl(-\hat\Lambda_0^r(t\,
#'     e^{\hat\alpha^\top z_i})\bigr), \quad
#'   S_{\text{de}}(t \mid z_i) = \exp\!\bigl(-\hat\Lambda_0^d(t\,
#'     e^{\hat\beta^\top z_i})\bigr).
#' }
#' The linear predictor \eqn{\hat\alpha^\top z_i} is a log time-acceleration
#' factor: the recurrence process for subject \eqn{i} runs on a time axis
#' scaled by \eqn{e^{\hat\alpha^\top z_i}} relative to the baseline.  Each
#' predictor contributes \eqn{\hat\alpha_j z_{ij}} to this log-scale factor,
#' so \eqn{e^{\hat\alpha_j z_{ij}}} is the multiplicative factor on the time
#' scale attributable to predictor \eqn{j} alone.  Values greater than 1
#' accelerate events (shorter times); values less than 1 decelerate them
#' (longer times).
#'
#' @param object An object of class \code{"swjm_cv"}.
#' @param newdata A numeric matrix or data frame of covariate values.
#'   Must have \code{p} columns named \code{x1}, \ldots, \code{xp} (if a
#'   data frame) or exactly \code{p} columns (if a matrix).
#' @param times Numeric vector of evaluation times (JFM only).  If
#'   \code{NULL}, the observed event times from the training data are used.
#' @param \dots Currently unused.
#'
#' @return An object of class \code{"swjm_pred"}, a list with:
#'   \describe{
#'     \item{S_re}{Matrix of readmission-free survival probabilities
#'       (rows = subjects, columns = \code{times}).}
#'     \item{S_de}{Matrix of death-free survival probabilities.}
#'     \item{times}{Numeric vector of evaluation times.}
#'     \item{lp_re}{Linear predictors for readmission
#'       (\eqn{\hat\alpha^\top z_i}).  For JSCM this is the log
#'       time-acceleration for the recurrence process.}
#'     \item{lp_de}{Linear predictors for death (\eqn{\hat\beta^\top z_i}).
#'       For JSCM this is the log time-acceleration for the terminal process.}
#'     \item{time_accel_re}{(JSCM only) \eqn{e^{\hat\alpha^\top z_i}}: the
#'       multiplicative factor by which the recurrence time axis is scaled
#'       relative to baseline.  \code{NULL} for JFM.}
#'     \item{time_accel_de}{(JSCM only) Analogous time-acceleration factor
#'       for the terminal process.  \code{NULL} for JFM.}
#'     \item{contrib_re}{Matrix of per-predictor contributions
#'       \eqn{\hat\alpha_j z_{ij}} (rows = subjects, columns = covariates).
#'       For JFM these are log-hazard contributions; for JSCM they are
#'       log time-acceleration contributions.}
#'     \item{contrib_de}{Analogous matrix for the terminal process.}
#'   }
#'
#' @examples
#' \donttest{
#' dat <- generate_data(n = 50, p = 5, scenario = 1, model = "jfm")
#' cv  <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
#'                     max_iter = 100)
#' newz <- matrix(rnorm(15), nrow = 3, ncol = 5)
#' pred <- predict(cv, newdata = newz)
#' plot(pred)
#'
#' dat_jscm <- generate_data(n = 50, p = 5, scenario = 1, model = "jscm")
#' cv_jscm  <- cv_stagewise(dat_jscm$data, model = "jscm", penalty = "coop",
#'                           max_iter = 500)
#' pred_jscm <- predict(cv_jscm, newdata = newz)
#' plot(pred_jscm)
#' }
#'
#' @export
predict.swjm_cv <- function(object, newdata, times = NULL, ...) {
  if (!inherits(object, "swjm_cv"))
    stop("'object' must be of class 'swjm_cv'")

  # Coerce newdata to a plain matrix of covariates
  if (is.data.frame(newdata)) {
    cov_cols <- paste0("x", seq_len(object$p))
    if (all(cov_cols %in% names(newdata))) {
      newdata <- as.matrix(newdata[, cov_cols])
    } else {
      newdata <- as.matrix(newdata)
    }
  } else {
    newdata <- as.matrix(newdata)
  }
  if (ncol(newdata) != object$p)
    stop(sprintf("'newdata' must have %d covariate column(s)", object$p))

  alpha  <- object$alpha
  beta   <- object$beta
  n_subj <- nrow(newdata)

  rn <- if (!is.null(rownames(newdata))) rownames(newdata) else
    paste0("Subject", seq_len(n_subj))
  xn <- paste0("x", seq_len(object$p))

  lp_re <- drop(newdata %*% alpha)
  lp_de <- drop(newdata %*% beta)
  names(lp_re) <- rn
  names(lp_de) <- rn

  contrib_re <- sweep(newdata, 2, alpha, `*`)
  contrib_de <- sweep(newdata, 2, beta,  `*`)
  rownames(contrib_re) <- rn;  colnames(contrib_re) <- xn
  rownames(contrib_de) <- rn;  colnames(contrib_de) <- xn

  if (object$model == "jfm") {
    if (is.null(object$baseline))
      stop("No baseline hazard stored; re-fit with cv_stagewise(model = 'jfm')")
    bl <- object$baseline
    if (is.null(times))
      times <- sort(unique(c(bl$tr, bl$td)))

    cum_r    <- c(0, cumsum(bl$lambda0_r))
    cum_d    <- c(0, cumsum(bl$lambda0_d))
    idx_r    <- findInterval(times, bl$tr) + 1L
    idx_d    <- findInterval(times, bl$td) + 1L
    Lambda_r <- cum_r[idx_r]
    Lambda_d <- cum_d[idx_d]

    S_re <- exp(outer(-exp(lp_re), Lambda_r))
    S_de <- exp(outer(-exp(lp_de), Lambda_d))
    cn <- paste0("t=", signif(times, 4))
    rownames(S_re) <- rn;  colnames(S_re) <- cn
    rownames(S_de) <- rn;  colnames(S_de) <- cn

    out <- list(
      S_re          = S_re,
      S_de          = S_de,
      times         = times,
      lp_re         = lp_re,
      lp_de         = lp_de,
      time_accel_re = NULL,
      time_accel_de = NULL,
      contrib_re    = contrib_re,
      contrib_de    = contrib_de,
      alpha         = alpha,
      beta          = beta,
      n_subj        = n_subj,
      p             = object$p,
      model         = object$model
    )
  } else {
    # JSCM: survival curves via Nelson-Aalen baseline on accelerated time scale.
    # For subject i at time t:
    #   S_re(t|z_i) = exp(-Λ_0^R(t · e^{α^T z_i}))
    #   S_de(t|z_i) = exp(-Λ_0^D(t · e^{β^T z_i}))
    bl <- object$baseline
    if (is.null(times))
      times <- sort(unique(c(bl$t_r, bl$t_d)))

    cum_r    <- c(0, cumsum(bl$lambda0_r))
    cum_d    <- c(0, cumsum(bl$lambda0_d))

    accel_re <- exp(lp_re)   # e^{α^T z_i}, length n_subj
    accel_de <- exp(lp_de)

    # For each subject i and each time t, evaluate Λ_0(t · accel_i)
    S_re <- matrix(NA_real_, nrow = n_subj, ncol = length(times))
    S_de <- matrix(NA_real_, nrow = n_subj, ncol = length(times))
    for (i in seq_len(n_subj)) {
      t_star_r <- times * accel_re[i]
      t_star_d <- times * accel_de[i]
      Lambda_r <- cum_r[findInterval(t_star_r, bl$t_r) + 1L]
      Lambda_d <- cum_d[findInterval(t_star_d, bl$t_d) + 1L]
      S_re[i, ] <- exp(-Lambda_r)
      S_de[i, ] <- exp(-Lambda_d)
    }
    cn <- paste0("t=", signif(times, 4))
    rownames(S_re) <- rn;  colnames(S_re) <- cn
    rownames(S_de) <- rn;  colnames(S_de) <- cn

    out <- list(
      S_re          = S_re,
      S_de          = S_de,
      times         = times,
      lp_re         = lp_re,
      lp_de         = lp_de,
      time_accel_re = accel_re,
      time_accel_de = accel_de,
      contrib_re    = contrib_re,
      contrib_de    = contrib_de,
      alpha         = alpha,
      beta          = beta,
      n_subj        = n_subj,
      p             = object$p,
      model         = object$model
    )
  }

  structure(out, class = "swjm_pred")
}


#' @export
print.swjm_pred <- function(x, ...) {
  cat("swjm predictions (", x$model, ")\n\n", sep = "")
  cat(sprintf("  %-24s %d\n", "Subjects:", x$n_subj))
  cat(sprintf("  %-24s %d\n", "Time points:", length(x$times)))
  cat(sprintf("  %-24s [%.4g, %.4g]\n", "Time range:",
              min(x$times), max(x$times)))
  if (x$model == "jscm") {
    cat("\n  Time-acceleration factors (exp(alpha^T z) for recurrence):\n")
    print(round(x$time_accel_re, 4))
    cat("\n  Time-acceleration factors (exp(beta^T z) for death):\n")
    print(round(x$time_accel_de, 4))
  }
  cat("\n  Use plot() to visualize survival curves and predictor contributions.\n")
  invisible(x)
}


#' Plot Predicted Survival Curves and Predictor Contributions
#'
#' Produces a figure for a \code{"swjm_pred"} object.  The layout depends on
#' the model:
#'
#' \strong{JFM} (four panels):
#' \enumerate{
#'   \item Readmission-free survival curves for all subjects, with the
#'     selected subject highlighted.
#'   \item Bar chart of readmission predictor contributions
#'     \eqn{\hat\alpha_j z_{ij}} (log-hazard scale).
#'   \item Mortality survival curves with the selected subject highlighted.
#'   \item Bar chart of death predictor contributions \eqn{\hat\beta_j z_{ij}}.
#' }
#'
#' \strong{JSCM} (four panels):
#' \enumerate{
#'   \item Recurrent-event survival curves (AFT model) with the selected
#'     subject highlighted.  The panel title shows the total acceleration
#'     factor \eqn{e^{\hat\alpha^\top z_i}}.
#'   \item Bar chart of recurrence log time-acceleration contributions
#'     \eqn{\hat\alpha_j z_{ij}}: positive = events sooner, negative = later.
#'   \item Mortality survival curves with the selected subject highlighted,
#'     total acceleration factor \eqn{e^{\hat\beta^\top z_i}} in the title.
#'   \item Bar chart of terminal-event log time-acceleration contributions
#'     \eqn{\hat\beta_j z_{ij}}.
#' }
#'
#' In all cases bars represent \emph{subject-specific} contributions
#' (coefficient \eqn{\times} covariate value), not bare coefficients, so the
#' display correctly reflects how much each predictor shifts the log-hazard
#' (JFM) or log time-acceleration (JSCM) for this particular subject.
#'
#' @param x An object of class \code{"swjm_pred"}.
#' @param which_subject Integer. Index of the subject to highlight (default
#'   \code{1}).
#' @param which_process Character. Which sub-model(s) to plot:
#'   \code{"both"} (default), \code{"readmission"}, or \code{"death"}.
#' @param threshold Non-negative numeric.  Only predictors whose absolute
#'   subject-specific contribution exceeds this value are shown in the bar
#'   chart.  The default \code{0} suppresses variables with exactly zero
#'   contribution (i.e., unselected variables whose coefficient is zero).
#'   Increase \code{threshold} to focus on the most impactful predictors.
#' @param \dots Currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' dat  <- generate_data(n = 50, p = 5, scenario = 1, model = "jfm")
#' cv   <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
#'                      max_iter = 100)
#' newz <- matrix(rnorm(15), nrow = 3, ncol = 5)
#' pred <- predict(cv, newdata = newz)
#' plot(pred, which_subject = 2)
#' plot(pred, which_subject = 2, threshold = 0.05)
#' }
#'
#' @export
plot.swjm_pred <- function(x,
                            which_subject  = 1L,
                            which_process  = c("both", "readmission", "death"),
                            threshold      = 0,
                            ...) {
  which_process <- match.arg(which_process)
  if (which_subject < 1L || which_subject > x$n_subj)
    stop(sprintf("'which_subject' must be between 1 and %d", x$n_subj))
  if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0)
    stop("'threshold' must be a single non-negative number")

  subj_label <- if (!is.null(rownames(x$contrib_re)))
    rownames(x$contrib_re)[which_subject] else paste0("Subject ", which_subject)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  both <- which_process == "both"

  # Both JFM and JSCM: 4 panels when both processes shown, 2 otherwise
  n_panels <- if (both) 4L else 2L
  if (n_panels == 4L)
    graphics::par(mfrow = c(2L, 2L), mar = c(4, 4, 3, 1))
  else
    graphics::par(mfrow = c(1L, 2L), mar = c(4, 4, 3, 1))

  # --- Survival panel ---
  .surv_panel <- function(S_mat, hl_col, title) {
    n <- nrow(S_mat)
    graphics::matplot(x$times, t(S_mat),
                      type = "l", lty = 1,
                      col  = grDevices::adjustcolor(seq_len(n), alpha.f = 0.25),
                      xlab = "Time", ylab = "S(t)",
                      main = title, ylim = c(0, 1))
    graphics::lines(x$times, S_mat[which_subject, ],
                    lwd = 2.5, col = hl_col)
    graphics::legend("topright",
                     legend = c("All subjects", subj_label),
                     col    = c("grey60", hl_col),
                     lwd    = c(1, 2.5), cex = 0.8, bty = "n")
  }

  # --- Contribution bar chart ---
  .contrib_panel <- function(contrib, main_title, ylab_expr) {
    cv   <- contrib[which_subject, ]
    keep <- abs(cv) > threshold
    if (!any(keep)) {
      graphics::plot.new()
      graphics::title(main = main_title)
      graphics::text(0.5, 0.5,
                     sprintf("No contributions exceed threshold (%g)", threshold),
                     adj = c(0.5, 0.5), col = "grey40")
      return(invisible(NULL))
    }
    cv  <- cv[keep]
    ord <- order(abs(cv), decreasing = TRUE)
    # Negative contribution = fewer/later events = better health -> blue (protective)
    # Positive contribution = more/sooner events = worse health -> red (harmful)
    bar_cols <- ifelse(cv[ord] >= 0, "tomato", "steelblue")
    rng  <- range(cv)
    pad  <- diff(rng) * 0.12
    if (pad == 0) pad <- abs(rng[1]) * 0.12 + 0.1
    ylim <- c(min(rng[1], 0) - pad, max(rng[2], 0) + pad)
    graphics::barplot(cv[ord],
                      names.arg = names(cv)[ord],
                      col       = bar_cols,
                      main      = main_title,
                      ylab      = ylab_expr,
                      ylim      = ylim,
                      las       = 2, cex.names = 0.75, cex.main = 0.85)
    graphics::abline(h = 0, lty = 2, col = "grey50")
  }

  if (x$model == "jfm") {
    surv_re_title  <- "Readmission-free survival"
    surv_de_title  <- "Survival"
    contrib_re_title <- paste0("Readmission log-hazard contributions\n", subj_label)
    contrib_de_title <- paste0("Death log-hazard contributions\n", subj_label)
    ylab_re <- expression(hat(alpha)[j] * z[j])
    ylab_de <- expression(hat(beta)[j] * z[j])
  } else {
    # JSCM: AFT interpretation — show time-acceleration factor in survival title
    accel_re     <- round(x$time_accel_re[which_subject], 3)
    accel_de     <- round(x$time_accel_de[which_subject], 3)
    dir_re       <- if (accel_re > 1) "events sooner" else "events later"
    dir_de       <- if (accel_de > 1) "death sooner"  else "death later"
    surv_re_title  <- paste0("Recurrent event survival (AFT)\n",
                              subj_label, "  [accel: ", accel_re, "x, ", dir_re, "]")
    surv_de_title  <- paste0("Survival (AFT)\n",
                              subj_label, "  [accel: ", accel_de, "x, ", dir_de, "]")
    contrib_re_title <- paste0("Recurrent event log time-acceleration\n", subj_label)
    contrib_de_title <- paste0("Terminal event log time-acceleration\n", subj_label)
    ylab_re <- expression(hat(alpha)[j] * z[j] ~~ (log ~ time-acceleration))
    ylab_de <- expression(hat(beta)[j] * z[j]  ~~ (log ~ time-acceleration))
  }

  if (which_process %in% c("both", "readmission")) {
    .surv_panel(x$S_re, "steelblue", surv_re_title)
    .contrib_panel(x$contrib_re, contrib_re_title, ylab_re)
  }
  if (which_process %in% c("both", "death")) {
    .surv_panel(x$S_de, "tomato", surv_de_title)
    .contrib_panel(x$contrib_de, contrib_de_title, ylab_de)
  }

  invisible(x)
}
