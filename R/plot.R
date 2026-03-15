# ============================================================================
# Plot methods for swjm_path and swjm_cv
# ============================================================================


#' Plot a Stagewise Coefficient Path
#'
#' Produces a glmnet-style plot of coefficient trajectories versus
#' \eqn{\log(\lambda)} along the stagewise regularization path.
#' Two panels are drawn by default: one for the readmission sub-model
#' (alpha) and one for the death sub-model (beta).  The top axis of
#' each panel shows the number of active variables at that lambda.
#'
#' @param x An object of class \code{"swjm_path"}.
#' @param log_lambda Logical. If \code{TRUE} (default), the x-axis is
#'   \eqn{\log(\lambda)}; otherwise raw \eqn{\lambda}.
#' @param which Character. Which sub-model(s) to plot:
#'   \code{"both"} (default), \code{"readmission"}, or \code{"death"}.
#' @param col Optional vector of colors, one per covariate. Recycled if
#'   shorter than \code{p}.
#' @param \dots Currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
#' fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
#'                      max_iter = 200)
#' plot(fit)
#' }
#'
#' @export
plot.swjm_path <- function(x,
                            log_lambda = TRUE,
                            which = c("both", "readmission", "death"),
                            col = NULL, ...) {
  which <- match.arg(which)
  p <- x$p

  dec_idx    <- extract_decreasing_indices(x$lambda)
  alpha_dec  <- x$alpha[, dec_idx, drop = FALSE]
  beta_dec   <- x$beta[, dec_idx, drop = FALSE]
  lambda_dec <- x$lambda[dec_idx]

  xv   <- if (log_lambda) log(lambda_dec) else lambda_dec
  xlab <- if (log_lambda) expression(log(lambda)) else expression(lambda)

  if (is.null(col))
    col <- seq_len(p)
  col <- rep_len(col, p)

  n_panels <- if (which == "both") 2L else 1L
  old_par  <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  # Extra top margin: line 0 = tick marks, line 1 = axis labels,
  # line 2.5 = "Active variables", line 4 = panel title.
  if (n_panels == 2L)
    graphics::par(mfrow = c(2L, 1L), mar = c(5, 4, 5, 2))
  else
    graphics::par(mar = c(5, 4, 5, 2))

  .draw_coef_panel <- function(coef_mat, process_label, active_counts) {
    ylim <- range(coef_mat, na.rm = TRUE)
    if (diff(ylim) == 0) ylim <- ylim + c(-1, 1) * 0.1
    xlim <- rev(range(xv))

    graphics::matplot(xv, t(coef_mat),
                      type = "l", lty = 1, lwd = 1.2, col = col,
                      xlab = xlab, ylab = "Coefficient",
                      main = "",          # title added separately below
                      xlim = xlim, ylim = ylim,
                      xaxt = "n")
    graphics::axis(1)
    graphics::abline(h = 0, lty = 2, col = "grey70")

    # Top axis: active variable count (line 0-1)
    at_vals <- pretty(xv, n = 5)
    at_vals <- at_vals[at_vals >= min(xv) & at_vals <= max(xv)]
    at_idx  <- vapply(at_vals, function(v) which.min(abs(xv - v)), integer(1))
    graphics::axis(3, at = at_vals, labels = active_counts[at_idx])
    # Axis label sits above the tick labels (line 2.5)
    graphics::mtext("Active variables", side = 3, line = 2.5, cex = 0.75)
    # Panel title sits above the axis label (line 4)
    graphics::title(main = paste0(process_label, "  [", x$model, "/",
                                  x$penalty, "]"),
                    line = 4, cex.main = 1)

    # Label nonzero variables at the right edge (smallest lambda)
    final_coefs <- coef_mat[, ncol(coef_mat)]
    nonzero_j   <- which(final_coefs != 0)
    if (length(nonzero_j) > 0) {
      graphics::text(min(xv), final_coefs[nonzero_j],
                     labels = paste0("x", nonzero_j),
                     cex = 0.65, adj = c(0, 0.5), col = col[nonzero_j])
    }
  }

  if (which %in% c("both", "readmission"))
    .draw_coef_panel(alpha_dec, "Readmission (alpha)",
                     colSums(alpha_dec != 0))
  if (which %in% c("both", "death"))
    .draw_coef_panel(beta_dec, "Death (beta)",
                     colSums(beta_dec != 0))

  invisible(x)
}


#' Plot Cross-Validation Results
#'
#' Plots the cross-validated estimating-equation score norms versus
#' \eqn{\log(\lambda)}, with separate lines for the readmission and death
#' components.  A vertical dashed line marks the lambda that minimizes the
#' combined norm.  The top axis shows the total number of active variables
#' along the path.
#'
#' @param x An object of class \code{"swjm_cv"}.
#' @param log_lambda Logical. If \code{TRUE} (default) the x-axis is
#'   \eqn{\log(\lambda)}; otherwise raw \eqn{\lambda}.
#' @param \dots Currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
#' cv  <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
#'                     max_iter = 200)
#' plot(cv)
#' }
#'
#' @export
plot.swjm_cv <- function(x, log_lambda = TRUE, ...) {
  xv   <- if (log_lambda) log(x$lambda_seq) else x$lambda_seq
  xlab <- if (log_lambda) expression(log(lambda)) else expression(lambda)

  ylim <- range(c(x$Scorenorm_crossfit,
                  x$Scorenorm_crossfit_re,
                  x$Scorenorm_crossfit_ce), na.rm = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  # Extra top margin: line 0 = ticks, line 1 = labels, line 2.5 = axis title,
  # line 4 = panel title.
  graphics::par(mar = c(5, 4, 5.5, 2))

  graphics::plot(xv, x$Scorenorm_crossfit,
                 type = "l", lwd = 2, col = "black",
                 xlab = xlab, ylab = "CV EE norm",
                 main = "",
                 ylim = ylim,
                 xlim = rev(range(xv)),
                 xaxt = "n")
  graphics::axis(1)
  graphics::lines(xv, x$Scorenorm_crossfit_re,
                  col = "steelblue", lwd = 1.5, lty = 2)
  graphics::lines(xv, x$Scorenorm_crossfit_ce,
                  col = "tomato",    lwd = 1.5, lty = 3)

  # Vertical line at optimal
  graphics::abline(v = xv[x$position_CF], lty = 2, col = "grey40")

  graphics::legend("topright",
                   legend = c("Combined", "Readmission", "Death"),
                   col    = c("black", "steelblue", "tomato"),
                   lty    = c(1, 2, 3),
                   lwd    = c(2, 1.5, 1.5),
                   cex    = 0.85, bty = "n")

  # Top axis: active variable count (line 0-1), label at line 2.5, title at 4
  at_vals <- pretty(xv, n = 5)
  at_vals <- at_vals[at_vals >= min(xv) & at_vals <= max(xv)]
  at_idx  <- vapply(at_vals, function(v) which.min(abs(xv - v)), integer(1))
  graphics::axis(3, at = at_vals, labels = x$n_active[at_idx])
  graphics::mtext("Active variables (total)", side = 3, line = 2.5, cex = 0.75)
  graphics::title(main = paste0("CV score norms  [", x$model, "/",
                                x$penalty, "]"),
                  line = 4.5, cex.main = 1)

  invisible(x)
}
