# ============================================================================
# Estimating equation functions for the Joint Scale-Change Model (JSCM)
# ============================================================================


#' Estimating Equation for Alpha (Recurrence) in the JSCM
#'
#' Evaluates the estimating equation for the recurrence coefficient vector
#' alpha in the joint scale-change model, using the C++ implementation
#' \code{am1}.
#'
#' @param Data2 A data frame in recurrent-event format.
#' @param a Numeric vector of current alpha coefficients.
#'
#' @return Numeric vector of score values, one per covariate.
#'
#' @keywords internal
jscm_ee_alpha <- function(Data2, a) {
  n <- length(unique(Data2$id))
  xi <- as.matrix(Data2[Data2$event == 0, 6:ncol(Data2)])
  rownames(xi) <- seq_len(nrow(xi))
  ti <- Data2$t.stop[Data2$event == 1]
  yi <- Data2$t.stop[Data2$event == 0]

  uids <- unique(Data2$id)
  Data_recur <- Data2[Data2$event == 1, ]
  list_recur <- vector("list", n)
  for (i in seq_len(n)) {
    list_recur[[i]] <- Data_recur$t.stop[Data_recur$id == uids[i]]
  }
  m <- as.integer(vapply(list_recur, length, integer(1)))
  w1 <- rep(1, length(yi))

  out <- am1(a, ti, yi, w1, xi, m)
  s <- if (is.matrix(out)) drop(out) else as.numeric(out)
  s
}


#' Estimating Equation for Beta (Death) in the JSCM
#'
#' Evaluates the estimating equation for the terminal event coefficient
#' vector beta in the joint scale-change model, using the C++ implementations
#' \code{temLog} and \code{reRate}.
#'
#' @param Data2 A data frame in recurrent-event format.
#' @param a Numeric vector of current alpha coefficients.
#' @param b Numeric vector of current beta coefficients.
#'
#' @return Numeric vector of score values, one per covariate.
#'
#' @keywords internal
jscm_ee_beta <- function(Data2, a, b) {
  n <- length(unique(Data2$id))
  xi <- as.matrix(Data2[Data2$event == 0, 6:ncol(Data2)])
  rownames(xi) <- seq_len(nrow(xi))
  ti <- Data2$t.stop[Data2$event == 1]
  yi <- Data2$t.stop[Data2$event == 0]

  uids <- unique(Data2$id)
  Data_recur <- Data2[Data2$event == 1, ]
  list_recur <- vector("list", n)
  for (i in seq_len(n)) {
    list_recur[[i]] <- Data_recur$t.stop[Data_recur$id == uids[i]]
  }
  m <- as.integer(vapply(list_recur, length, integer(1)))
  w1 <- rep(1, length(yi))
  di <- Data2$status[Data2$event == 0]

  # estimate of frailty term
  texa <- log(ti) + as.matrix(Data2[Data2$event == 1, -(1:5)]) %*% a
  yexa <- log(yi) + xi %*% a
  rate <- c(reRate(texa, rep(yexa, m), rep(w1, m), yexa))
  Lam <- exp(-rate)
  R <- m / Lam
  numAdj <- 1e-3
  if (numAdj > min(Lam)) {
    numAdj <- numAdj * min(Lam)
  }
  R2 <- (m + numAdj) / (Lam + numAdj)

  out <- temLog(b, b, xi, yi, R2, di, w1)
  s <- if (is.matrix(out)) drop(out) else as.numeric(out)
  s
}
