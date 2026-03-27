# ============================================================================
# Shared utility functions used by both JFM and JSCM implementations
# ============================================================================

#' Interpolate Coefficient Vectors Along a Lambda Sequence
#'
#' Given a decreasing lambda sequence from a fitted path and a new set of lambda
#' values, linearly interpolates the corresponding coefficient vectors.
#'
#' @param lambda Numeric vector of decreasing lambda values from the fitted path.
#' @param s Numeric vector of new lambda values at which to interpolate.
#'
#' @return A list with components:
#'   \item{left}{Integer vector of left neighbor indices.}
#'   \item{right}{Integer vector of right neighbor indices.}
#'   \item{frac}{Numeric vector of interpolation fractions. The interpolated
#'     value is \code{frac * path[left] + (1 - frac) * path[right]}.}
#'
#' @keywords internal
lambda_interp <- function(lambda, s) {
  if (length(lambda) == 1L) {
    nums <- length(s)
    left <- rep(1L, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    k <- length(lambda)
    sfrac <- (lambda[1] - s) / (lambda[1] - lambda[k])
    lambda_norm <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    sfrac[sfrac < min(lambda_norm)] <- min(lambda_norm)
    sfrac[sfrac > max(lambda_norm)] <- max(lambda_norm)
    coord <- approx(lambda_norm, seq_along(lambda_norm), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda_norm[right]) / (lambda_norm[left] - lambda_norm[right])
    sfrac[left == right] <- 1
    sfrac[abs(lambda_norm[left] - lambda_norm[right]) < .Machine$double.eps] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}


#' Compute the Cooperative Lasso Norm
#'
#' For each pair (theta_j, theta_(j+p)), computes the cooperative lasso norm:
#' L2 norm if signs agree, L1 norm if signs disagree.
#'
#' @param theta Numeric vector of length 2p (concatenation of alpha and beta).
#' @param p Integer, number of covariates.
#'
#' @return Scalar cooperative lasso norm value.
#'
#' @keywords internal
coop_norm <- function(theta, p) {
  h <- numeric(p)
  for (i in seq_len(p)) {
    pair <- c(theta[i], theta[p + i])
    if (pair[1] * pair[2] >= 0) {
      h[i] <- sqrt(sum(pair^2))
    } else {
      h[i] <- sum(abs(pair))
    }
  }
  sum(h)
}


#' Count Leading Zeros After Decimal Point
#'
#' Used for adaptive step size computation. Returns the position of the first
#' non-zero digit after the decimal point.
#'
#' @param num A numeric scalar.
#'
#' @return Integer count of leading zeros after the decimal point.
#'
#' @keywords internal
count_digits <- function(num) {
  num_str <- as.character(num)
  decimal_pos <- gregexpr("\\.", num_str)[[1]][1]
  if (decimal_pos == -1) return(0L)
  after_decimal <- substring(num_str, decimal_pos + 1)
  first_nonzero_pos <- gregexpr("[1-9]", after_decimal)[[1]][1]
  if (first_nonzero_pos == -1) return(0L)
  as.integer(first_nonzero_pos)
}


#' Create Stratified K-Fold Splits
#'
#' Randomly assigns a vector of IDs into K approximately equal-sized folds.
#'
#' @param ids Vector of subject IDs.
#' @param K Integer number of folds.
#'
#' @return A list of length K, where each element contains the IDs assigned
#'   to that fold.
#'
#' @keywords internal
create_folds <- function(ids, K) {
  n <- length(ids)
  folds <- sample(rep(seq_len(K), length.out = n))
  split(ids, folds)
}


#' Extract Common Data Components from a Recurrent-Event Data Frame
#'
#' Parses a standard recurrent-event data frame into the components needed
#' by the JFM estimating equations: covariate lists, event times, at-risk
#' indicators, etc.
#'
#' @param Data2 A data frame in the standard recurrent-event format.
#'
#' @return A list with components:
#'   \item{Z}{List of covariate matrices, one per subject.}
#'   \item{n}{Number of unique subjects.}
#'   \item{p}{Number of covariates.}
#'   \item{td}{Death times.}
#'   \item{td.id}{Subject IDs for death times.}
#'   \item{d_td}{Table of death time frequencies.}
#'   \item{tr}{Recurrent event times.}
#'   \item{tr.id}{Subject IDs for recurrent events.}
#'   \item{d_tr}{Table of recurrent event frequencies.}
#'   \item{Y}{Composite censoring/death times (one per subject).}
#'   \item{STATUS}{Death indicator at composite time (one per subject).}
#'   \item{list_recur}{List of recurrent event times per subject.}
#'   \item{num_recur}{Integer vector of recurrent event counts per subject.}
#'   \item{t.start}{All interval start times.}
#'   \item{I}{All subject IDs (matching rows of the data frame).}
#'
#' @keywords internal
extract_data_components <- function(Data2) {
  p <- ncol(Data2) - 5L
  uids <- unique(Data2$id)
  n <- length(uids)

  # Covariate matrices per subject — vectorized via split()
  cov_cols <- as.matrix(Data2[, 6:ncol(Data2), drop = FALSE])
  row_groups <- split(seq_len(nrow(Data2)), Data2$id)
  Z <- lapply(row_groups[as.character(uids)], function(rows) cov_cols[rows, , drop = FALSE])

  # Death events
  is_death <- Data2$status == 1
  td <- Data2$t.stop[is_death]
  td.id <- Data2$id[is_death]
  d_td <- table(td)

  # Recurrent events
  is_recur <- Data2$event == 1
  tr <- Data2$t.stop[is_recur]
  tr.id <- Data2$id[is_recur]
  d_tr <- table(tr)

  # Composite censoring/death times and status (one per subject, event==0 rows)
  is_last <- Data2$event == 0
  Y <- Data2$t.stop[is_last]
  STATUS <- Data2$status[is_last]

  # Recurrent events per subject — vectorized via split()
  recur_stops <- Data2$t.stop[is_recur]
  recur_ids   <- Data2$id[is_recur]
  list_recur_raw <- split(recur_stops, recur_ids)
  list_recur <- vector("list", n)
  for (i in seq_len(n)) {
    uid_char <- as.character(uids[i])
    list_recur[[i]] <- if (uid_char %in% names(list_recur_raw)) list_recur_raw[[uid_char]] else numeric(0)
  }
  num_recur <- vapply(list_recur, length, integer(1))

  t.start <- Data2$t.start
  I <- Data2$id

  list(
    Z = Z, n = n, p = p,
    td = td, td.id = td.id, d_td = d_td,
    tr = tr, tr.id = tr.id, d_tr = d_tr,
    Y = Y, STATUS = STATUS,
    list_recur = list_recur, num_recur = num_recur,
    t.start = t.start, I = I
  )
}


#' Extract Decreasing Lambda Path
#'
#' Post-processes a raw lambda sequence from the stagewise algorithm to keep
#' only strictly decreasing values (via running minimum and deduplication).
#'
#' @param lambda Numeric vector of raw lambda values from stagewise iterations.
#' @param tol_digits Integer, number of digits for rounding in deduplication.
#'
#' @return An integer vector of indices into the original lambda vector
#'   corresponding to the strictly decreasing subsequence.
#'
#' @keywords internal
extract_decreasing_indices <- function(lambda, tol_digits = 6L) {
  m <- cummin(lambda)
  m_round <- round(m, tol_digits)
  which(!duplicated(m_round))
}
