# ============================================================================
# Shared utility functions used by both JFM and JSCM implementations
# ============================================================================


#' Prepare Data for swjm Functions
#'
#' Expands factor and character covariate columns into dummy (indicator)
#' variables. The first 5 columns (id, t.start, t.stop, event, status) are
#' left unchanged. For each factor/character column with L levels, L-1 dummy
#' columns are created (reference level dropped). A message is printed when
#' expansion occurs.
#'
#' @param data A data frame.
#' @param caller Character string naming the calling function (for messages).
#'
#' @return A data frame with factor/character covariates replaced by numeric
#'   dummy columns.
#'
#' @keywords internal
prepare_data <- function(data, caller = "swjm") {
  if (!is.data.frame(data) || ncol(data) <= 5L) return(data)

  meta_cols <- data[, 1:5, drop = FALSE]
  cov_df    <- data[, 6:ncol(data), drop = FALSE]

  needs_expand <- vapply(cov_df, function(x) {
    is.factor(x) || is.character(x) || is.logical(x)
  }, logical(1))

  if (!any(needs_expand)) return(data)

  expanded_names <- names(cov_df)[needs_expand]
  message(caller, ": expanding non-numeric covariates to dummy variables: ",
          paste(expanded_names, collapse = ", "))

  # Convert characters to factors
  for (j in which(needs_expand)) {
    if (is.character(cov_df[[j]]) || is.logical(cov_df[[j]])) {
      cov_df[[j]] <- factor(cov_df[[j]])
    }
  }
  
  # For columns that were already numeric, model.matrix keeps them as-is.
  # For factors, it creates L-1 dummy columns (treatment contrasts).
  # We need to remove the original factor columns' single column and
  # replace with the dummies.

  # Simpler approach: build model matrix from scratch with proper contrasts
  # Use contrasts that drop the first level (default treatment contrasts)
  frm <- reformulate(names(cov_df))
  mm <- model.matrix(frm, data = cov_df)
  # Remove intercept column
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

  result <- cbind(meta_cols, as.data.frame(mm))
  rownames(result) <- rownames(data)
  result
}


#' Validate Recurrent-Event Data Frame
#'
#' Checks that a data frame has the required structure for swjm functions:
#' columns id, t.start, t.stop, event, status, and at least one covariate.
#' Produces informative error messages for each violated requirement.
#'
#' @param data Object to validate.
#' @param caller Character string naming the calling function (for error messages).
#'
#' @return Invisibly returns \code{TRUE} if valid; throws an error otherwise.
#'
#' @keywords internal
validate_data <- function(data, caller = "swjm") {
  problems <- character(0)

  if (!is.data.frame(data)) {
    stop(caller, ": 'data' must be a data frame, got ", class(data)[1], ".",
         call. = FALSE)
  }

  required <- c("id", "t.start", "t.stop", "event", "status")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    problems <- c(problems,
      paste0("Missing required columns: ", paste(missing_cols, collapse = ", "),
             ". Expected: id, t.start, t.stop, event, status, x1, ..., xp."))
  }

  p <- ncol(data) - 5L
  if (p < 1) {
    problems <- c(problems,
      paste0("Data must have at least 6 columns (5 required + covariates), ",
             "but has ", ncol(data), "."))
  }

  if (length(problems) == 0 && all(required %in% names(data))) {
    # Check column order: first 5 must be id, t.start, t.stop, event, status
    if (!identical(names(data)[1:5], required)) {
      problems <- c(problems,
        paste0("First 5 columns must be id, t.start, t.stop, event, status (in order), ",
               "but got: ", paste(names(data)[1:5], collapse = ", "), "."))
    }

    # Check numeric types
    for (col in c("t.start", "t.stop")) {
      if (!is.numeric(data[[col]])) {
        problems <- c(problems,
          paste0("Column '", col, "' must be numeric, got ", class(data[[col]])[1], "."))
      }
    }

    # Check event and status are 0/1
    for (col in c("event", "status")) {
      vals <- unique(data[[col]])
      if (!all(vals %in% c(0, 1))) {
        problems <- c(problems,
          paste0("Column '", col, "' must contain only 0 and 1, ",
                 "got values: ", paste(sort(unique(vals)), collapse = ", "), "."))
      }
    }

    # Check t.start <= t.stop
    if (is.numeric(data$t.start) && is.numeric(data$t.stop)) {
      bad <- sum(data$t.start > data$t.stop)
      if (bad > 0) {
        problems <- c(problems,
          paste0(bad, " row(s) have t.start > t.stop."))
      }
    }

    # Check each subject has exactly one terminal row (event == 0)
    n_terminal <- tapply(data$event == 0, data$id, sum)
    bad_subj <- names(n_terminal)[n_terminal != 1]
    if (length(bad_subj) > 0) {
      n_show <- min(length(bad_subj), 5)
      problems <- c(problems,
        paste0(length(bad_subj), " subject(s) do not have exactly one terminal row ",
               "(event == 0): ", paste(bad_subj[1:n_show], collapse = ", "),
               if (length(bad_subj) > n_show) ", ..." else "", "."))
    }

    # Check covariates are numeric
    cov_cols <- 6:ncol(data)
    non_numeric <- names(data)[cov_cols][!vapply(data[cov_cols], is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
      problems <- c(problems,
        paste0("Covariate columns must be numeric. Non-numeric: ",
               paste(non_numeric, collapse = ", "), "."))
    }

    # Check for NA/NaN/Inf
    has_na <- vapply(data, function(x) any(is.na(x) | is.nan(x)), logical(1))
    if (any(has_na)) {
      problems <- c(problems,
        paste0("Data contains NA/NaN values in columns: ",
               paste(names(data)[has_na], collapse = ", "), "."))
    }
    has_inf <- vapply(data[cov_cols], function(x) any(is.infinite(x)), logical(1))
    if (any(has_inf)) {
      problems <- c(problems,
        paste0("Data contains Inf values in covariate columns: ",
               paste(names(data)[cov_cols][has_inf], collapse = ", "), "."))
    }
  }

  if (length(problems) > 0) {
    stop(caller, ": invalid data frame.\n  ",
         paste(problems, collapse = "\n  "),
         call. = FALSE)
  }

  invisible(TRUE)
}

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
