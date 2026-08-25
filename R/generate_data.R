#' Generate Simulated Data for Joint Models
#'
#' Unified interface that dispatches to model-specific data generation
#' functions for the joint frailty model (JFM) or joint scale-change
#' model (JSCM).
#'
#' @param n Integer. Number of subjects.
#' @param p Integer. Number of covariates (should be a multiple of 10 for
#'   scenarios 1--3).
#' @param scenario Integer. Scenario for true coefficient configuration
#'   (1, 2, 3, or other for a simple default).
#' @param model Character. Either \code{"jfm"} for the joint frailty model
#'   or \code{"jscm"} for the joint scale-change model.
#' @param ... Additional arguments passed to the model-specific function.
#'   For JFM: \code{b}, \code{lambda0_d}, \code{lambda0_r}, \code{gamma_frailty},
#'   \code{cov_type}.
#'   For JSCM: \code{b}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{data}{A data frame in recurrent-event format with columns
#'       \code{id}, \code{t.start}, \code{t.stop}, \code{event},
#'       \code{status}, and covariate columns \code{x1}, \ldots, \code{xp}.}
#'     \item{alpha_true}{Numeric vector of true alpha coefficients.}
#'     \item{beta_true}{Numeric vector of true beta coefficients.}
#'   }
#'
#' @examples
#' # JFM data with 30 subjects and 10 covariates
#' dat_jfm <- generate_data(n = 30, p = 10, scenario = 1, model = "jfm")
#' head(dat_jfm$data)
#'
#' # JSCM data with 30 subjects and 10 covariates
#' dat_jscm <- generate_data(n = 30, p = 10, scenario = 1, model = "jscm")
#' head(dat_jscm$data)
#'
#' @export
generate_data <- function(n, p, scenario = 1, model = c("jfm", "jscm"), ...) {
  model <- match.arg(model)
  if (model == "jfm") {
    generate_data_jfm(n = n, p = p, scenario = scenario, ...)
  } else {
    generate_data_jscm(n = n, p = p, scenario = scenario, ...)
  }
}


#' Generate Simulated Data for the Joint Frailty Model (JFM)
#'
#' Generates recurrent-event and terminal-event data under a Cox-type
#' joint frailty model. Ported from
#' \code{Data_Generation_time_dependent_new()}.
#'
#' Internally the simulation uses the Rondeau et al. (2007) convention
#' where \code{alpha} governs death and \code{beta} governs recurrence.
#' The returned \code{alpha_true} and \code{beta_true} are relabelled to
#' match the package-wide convention:
#' \itemize{
#'   \item \code{alpha_true}: recurrence (readmission) coefficients.
#'   \item \code{beta_true}: terminal (death) coefficients.
#' }
#'
#' Within each subject the covariates are regenerated at each gap time,
#' yielding time-dependent covariates. Censoring times are
#' \code{Uniform(1, b)}.
#'
#' @param n Integer. Number of subjects.
#' @param p Integer. Number of covariates.
#' @param scenario Integer. Scenario (1, 2, 3, or other).
#' @param b Numeric. Upper bound of the censoring uniform distribution
#'   (default 6.50).
#' @param lambda0_d Numeric. Baseline hazard rate for the terminal event
#'   (default 0.041).
#' @param lambda0_r Numeric. Baseline hazard rate for recurrent events
#'   (default 1).
#' @param gamma_frailty Numeric. Frailty variance parameter. When positive,
#'   a subject-specific frailty \eqn{Z_i \sim \text{Gamma}(1/\gamma, 1/\gamma)}
#'   is drawn for each subject and multiplies both hazard rates. When 0
#'   (default), no frailty is used (\eqn{Z_i = 1}).
#' @param cov_type Character. How time-varying covariates are generated:
#'   \code{"internal"} (default) redraws covariates at each recurrent event;
#'   \code{"external"} changes covariates at predetermined Poisson times
#'   independent of the event process (Kalbfleisch-compatible);
#'   \code{"fixed"} draws one covariate vector per subject that never changes.
#'
#' @return A list with components:
#'   \describe{
#'     \item{data}{Data frame with columns \code{id}, \code{t.start},
#'       \code{t.stop}, \code{event}, \code{status}, \code{x1}, \ldots,
#'       \code{xp}.}
#'     \item{alpha_true}{True alpha (terminal) coefficients.}
#'     \item{beta_true}{True beta (recurrence) coefficients.}
#'   }
#'
#' @examples
#' dat <- generate_data_jfm(n = 30, p = 10, scenario = 1)
#' head(dat$data)
#' dat$alpha_true
#' dat$beta_true
#'
#' @export
generate_data_jfm <- function(n, p, scenario = 1, b = 6.50,
                               lambda0_d = 0.041, lambda0_r = 1,
                               gamma_frailty = 0,
                               cov_type = c("internal", "external", "fixed")) {
  cov_type <- match.arg(cov_type)
  S <- scenario
  theta <- 1
  alpha.star <- numeric(p)
  beta.star  <- numeric(p)

  # Helper: indices for the k-th 10% block (1-indexed), safe for any p >= 1
  .blk <- function(k) {
    lo <- floor(p * (k - 1) / 10) + 1L
    hi <- floor(p * k / 10)
    if (hi >= lo) lo:hi else integer(0)
  }

  if (S == 1) {
    beta.star[.blk(1)]  <-  1.1;  alpha.star[.blk(1)]  <-  0.1
    beta.star[.blk(2)]  <- -1.1;  alpha.star[.blk(2)]  <- -0.1
    beta.star[.blk(3)]  <-  0.1;  alpha.star[.blk(3)]  <-  1.1
    beta.star[.blk(4)]  <- -0.1;  alpha.star[.blk(4)]  <- -1.1
    beta.star[.blk(9)]  <-  1;    alpha.star[.blk(9)]  <-  1
    beta.star[.blk(10)] <- -1;    alpha.star[.blk(10)] <- -1

  } else if (S == 2) {
    beta.star[.blk(1)]  <-  1;    alpha.star[.blk(1)]  <-  1
    beta.star[.blk(2)]  <- -1;    alpha.star[.blk(2)]  <- -1
    beta.star[.blk(3)]  <- -1;    alpha.star[.blk(3)]  <-  1
    beta.star[.blk(4)]  <-  1;    alpha.star[.blk(4)]  <- -1
    beta.star[.blk(9)]  <-  1;    alpha.star[.blk(10)] <-  1

  } else if (S == 3) {
    beta.star[.blk(1)]  <-  1;    alpha.star[.blk(1)]  <- -1
    beta.star[.blk(2)]  <- -1;    alpha.star[.blk(2)]  <-  1
    beta.star[.blk(3)]  <-  1.5;  alpha.star[.blk(3)]  <- -1
    beta.star[.blk(4)]  <- -1;    alpha.star[.blk(4)]  <-  1.5
    beta.star[.blk(9)]  <-  1
    alpha.star[.blk(10)] <- 1
  } else {
    beta.star <- alpha.star <- c(-0.5, 0.5)
  }

  alpha <- alpha.star
  beta  <- beta.star

  # Pre-allocate list of per-subject results (avoids quadratic vector growing)
  subj_list <- vector("list", n)

  censors <- runif(n, 1, b)
  frailties <- if (gamma_frailty > 0) {
    rgamma(n, shape = 1 / gamma_frailty, rate = 1 / gamma_frailty)
  } else {
    rep(1, n)
  }

  for (j in 1:n) {
    censor <- censors[j]
    Z_i    <- frailties[j]

    # Generate covariate schedule based on cov_type
    if (cov_type == "fixed") {
      # Single covariate drawn once per subject
      z_fixed <- rnorm(p)
      cov_break_times <- c(0, censor)
      cov_values <- list(z_fixed)
    } else if (cov_type == "external") {
      # Covariates change at predetermined Poisson times (independent of events)
      n_changes <- rpois(1, lambda = 3)
      ch_times <- if (n_changes > 0) sort(runif(n_changes, 0, censor)) else numeric(0)
      cov_break_times <- c(0, ch_times, censor)
      cov_values <- lapply(seq_len(length(cov_break_times) - 1), function(i) rnorm(p))
    } else {
      # "internal": covariates redrawn at each recurrent event (legacy behavior)
      cov_break_times <- NULL  # not used; z drawn inside the event loop
      cov_values <- NULL
    }

    Xk    <- 0
    delta <- 0
    max_events <- 200L
    X_buf   <- numeric(max_events)
    ev_buf  <- integer(max_events)
    st_buf  <- integer(max_events)
    cov_buf <- matrix(NA_real_, max_events, p)
    k <- 0L

    if (cov_type == "internal") {
      # Legacy: redraw covariates at each event
      while (Xk < censor && delta == 0L) {
        z <- rnorm(p)
        T_gap <- rexp(1, rate = Z_i * lambda0_r * exp(drop(beta %*% z)))
        D_gap <- rexp(1, rate = Z_i * lambda0_d * exp(drop(alpha %*% z)))

        Xk <- Xk + min(T_gap, D_gap)
        if (D_gap < T_gap && Xk < censor) delta <- 1L else delta <- 0L

        k <- k + 1L
        if (k > max_events) {
          max_events <- max_events * 2L
          X_buf   <- c(X_buf, numeric(max_events - length(X_buf)))
          ev_buf  <- c(ev_buf, integer(max_events - length(ev_buf)))
          st_buf  <- c(st_buf, integer(max_events - length(st_buf)))
          cov_buf <- rbind(cov_buf, matrix(NA_real_, max_events - nrow(cov_buf), p))
        }
        X_buf[k]     <- Xk
        ev_buf[k]    <- 1L - delta
        st_buf[k]    <- delta
        cov_buf[k, ] <- z
      }
    } else {
      # "fixed" or "external": process covariate intervals
      n_intv <- length(cov_break_times) - 1L
      for (ci in seq_len(n_intv)) {
        z <- cov_values[[ci]]
        intv_end <- cov_break_times[ci + 1L]
        intv_start <- if (k > 0L) X_buf[k] else 0

        while (Xk < intv_end && delta == 0L) {
          T_gap <- rexp(1, rate = Z_i * lambda0_r * exp(drop(beta %*% z)))
          D_gap <- rexp(1, rate = Z_i * lambda0_d * exp(drop(alpha %*% z)))
          gap <- min(T_gap, D_gap)

          if (Xk + gap >= intv_end) break  # no event before interval ends

          Xk <- Xk + gap
          if (Xk >= censor) break
          if (D_gap < T_gap) delta <- 1L else delta <- 0L

          k <- k + 1L
          if (k > max_events) {
            max_events <- max_events * 2L
            X_buf   <- c(X_buf, numeric(max_events - length(X_buf)))
            ev_buf  <- c(ev_buf, integer(max_events - length(ev_buf)))
            st_buf  <- c(st_buf, integer(max_events - length(st_buf)))
            cov_buf <- rbind(cov_buf, matrix(NA_real_, max_events - nrow(cov_buf), p))
          }
          X_buf[k]     <- Xk
          ev_buf[k]    <- 1L - delta
          st_buf[k]    <- delta
          cov_buf[k, ] <- z
          intv_start <- Xk
        }
        if (delta == 1L) break
        Xk <- intv_end  # advance to next interval
      }
    }

    if (k == 0L) {
      # No events: single censoring row
      k <- 1L
      X_buf[1]     <- censor
      ev_buf[1]    <- 0L
      st_buf[1]    <- 0L
      cov_buf[1, ] <- if (cov_type == "internal") rnorm(p) else cov_values[[1]]
    }

    X_buf  <- X_buf[1:k]
    ev_buf <- ev_buf[1:k]
    st_buf <- st_buf[1:k]
    cov_buf <- cov_buf[1:k, , drop = FALSE]

    if (delta != 1L) {
      X_buf[k]  <- censor
      ev_buf[k] <- 0L
    }

    sub_t_start <- c(0, X_buf[-k])
    subj_list[[j]] <- list(
      id      = rep(j, k),
      t.start = sub_t_start,
      t.stop  = X_buf,
      event   = ev_buf,
      status  = st_buf,
      cov     = cov_buf
    )
  }

  # Collapse all subjects at once
  subj_list <- subj_list[!vapply(subj_list, is.null, logical(1))]
  id         <- unlist(lapply(subj_list, `[[`, "id"))
  t.start    <- unlist(lapply(subj_list, `[[`, "t.start"))
  t.stop     <- unlist(lapply(subj_list, `[[`, "t.stop"))
  event      <- unlist(lapply(subj_list, `[[`, "event"))
  status     <- unlist(lapply(subj_list, `[[`, "status"))
  cov.matrix <- do.call(rbind, lapply(subj_list, `[[`, "cov"))

  colnames(cov.matrix) <- paste0("x", 1:p)
  Data2 <- data.frame(
    id = id, t.start = t.start, t.stop = t.stop,
    event = event, status = status, cov.matrix
  )
  rownames(Data2) <- 1:nrow(Data2)

  # Relabel to package convention: alpha = recurrence, beta = death
  # (internally alpha.star drove death and beta.star drove recurrence)
  list(data = Data2, alpha_true = beta.star, beta_true = alpha.star,
       gamma_frailty = gamma_frailty)
}


#' Generate Simulated Data for the Joint Scale-Change Model (JSCM)
#'
#' Generates recurrent-event and terminal-event data under an AFT-type
#' joint scale-change model using \code{\link[reReg]{simGSC}}.
#' Ported from \code{Data_gen_reReg()}.
#'
#' In the JSCM convention:
#' \itemize{
#'   \item \code{alpha} governs the recurrence process.
#'   \item \code{beta} governs the terminal (death) process.
#' }
#'
#' Covariates are drawn from \code{Uniform(-1, 1)}. A gamma frailty with
#' shape = 4, scale = 1/4 is used. Censoring times are
#' \code{Uniform(0, b)}.
#'
#' @param n Integer. Number of subjects.
#' @param p Integer. Number of covariates.
#' @param scenario Integer. Scenario (1, 2, 3, or other).
#' @param b Numeric. Upper bound of the censoring uniform distribution
#'   (default 4).
#'
#' @return A list with components:
#'   \describe{
#'     \item{data}{Object returned by \code{\link[reReg]{simGSC}} (a
#'       data frame with recurrent-event structure).}
#'     \item{alpha_true}{True alpha (recurrence) coefficients.}
#'     \item{beta_true}{True beta (terminal) coefficients.}
#'   }
#'
#' @examples
#' \donttest{
#' dat <- generate_data_jscm(n = 30, p = 10, scenario = 1)
#' head(dat$data)
#' dat$alpha_true
#' dat$beta_true
#' }
#'
#' @importFrom reReg simGSC
#' @export
generate_data_jscm <- function(n, p, scenario = 1, b = 4) {
  S <- scenario

  # The scenarios below place their nonzero coefficients at fixed positions
  # within the first ten covariates, so p must leave room for them. With a
  # smaller p, assigning position 10 would silently extend the coefficient
  # vector past the covariate matrix, and the mismatch would surface much
  # later as an opaque message from reReg::simGSC about alpha not matching
  # the number of covariates.
  if (p < 10) {
    stop("'p' must be at least 10 for the joint scale-change scenarios, ",
         "which place their nonzero coefficients within the first ten ",
         "covariates (got p = ", p, ").", call. = FALSE)
  }

  alpha.star <- numeric(p)
  beta.star  <- numeric(p)

  if (S == 1) {
    # first 10% large alpha and small beta (+)
    alpha.star[1:(10 * 0.1)] <- 1.1
    beta.star[1:(10 * 0.1)] <- 0.1
    # 10%-20% large alpha and small beta (-)
    alpha.star[(10 * 0.1 + 1):(10 * 0.2)] <- -1.1
    beta.star[(10 * 0.1 + 1):(10 * 0.2)] <- -0.1
    # 20-30% small alpha and large beta (+)
    alpha.star[(10 * 0.2 + 1):(10 * 0.3)] <- 0.1
    beta.star[(10 * 0.2 + 1):(10 * 0.3)] <- 1.1
    # 30%-40% small alpha and large beta (-)
    alpha.star[(10 * 0.3 + 1):(10 * 0.4)] <- -0.1
    beta.star[(10 * 0.3 + 1):(10 * 0.4)] <- -1.1
    # 80-90% similar magnitude (+)
    alpha.star[(10 * 0.8 + 1):(10 * 0.9)] <- 1
    beta.star[(10 * 0.8 + 1):(10 * 0.9)] <- 1
    # 90-100% similar magnitude (-)
    alpha.star[(10 * 0.9 + 1):(10 * 1)] <- -1
    beta.star[(10 * 0.9 + 1):(10 * 1)] <- -1

  } else if (S == 2) {
    # give the 10% to 20% both strong signal and 20% to 30% both weak signal
    alpha.star[(10 * 0.1 + 1):(10 * 0.2)] <- -1
    alpha.star[(10 * 0.2 + 1):(10 * 0.3)] <- -1

    # first 10%
    alpha.star[1:(10 * 0.1)] <- 1

    # 30% to 40%
    alpha.star[(10 * 0.3 + 1):(10 * 0.4)] <- 1

    # 80% to 90%
    alpha.star[(10 * 0.8 + 1):(10 * 0.9)] <- 1

    # 0-10%
    beta.star[1:(10 * 0.1)] <- 1
    # 10%-20%
    beta.star[(10 * 0.1 + 1):(10 * 0.2)] <- -1
    # 20%-30%
    beta.star[(10 * 0.2 + 1):(10 * 0.3)] <- 1
    # 30%-40%
    beta.star[(10 * 0.3 + 1):(10 * 0.4)] <- -1
    # 90% to 100%
    beta.star[(10 * 0.9 + 1):10] <- 1

  } else if (S == 3) {
    # 10%: alpha:1 beta: -1
    alpha.star[1:(10 * 0.1)] <- 1
    beta.star[1:(10 * 0.1)] <- -1
    # 10-20%: alpha:-1 beta: 1
    alpha.star[(10 * 0.1 + 1):(10 * 0.2)] <- -1
    beta.star[(10 * 0.1 + 1):(10 * 0.2)] <- 1
    # 20-30%: alpha:2 beta: -1
    alpha.star[(10 * 0.2 + 1):(10 * 0.3)] <- 1.5
    beta.star[(10 * 0.2 + 1):(10 * 0.3)] <- -1
    # 30-40%: alpha:-1 beta: 2
    alpha.star[(10 * 0.3 + 1):(10 * 0.4)] <- -1
    beta.star[(10 * 0.3 + 1):(10 * 0.4)] <- 1.5
    # 80-90%: alpha:1 beta: 0
    alpha.star[(10 * 0.8 + 1):(10 * 0.9)] <- 1
    beta.star[(10 * 0.8 + 1):(10 * 0.9)] <- 0
    # 90-100%: alpha:0 beta: 1
    alpha.star[(10 * 0.9 + 1):10] <- 0
    beta.star[(10 * 0.9 + 1):10] <- 1
  } else {
    alpha.star <- beta.star <- c(-0.5, 0.5)
  }

  alpha <- beta <- alpha.star
  eta   <- theta <- beta.star

  X <- matrix(runif(n * p, -1, 1), n, p)
  colnames(X) <- paste0("x", 1:p)

  C <- runif(n, min = 0, max = b)

  para <- list(alpha = alpha, beta = beta, eta = eta, theta = theta)

  gamma <- rgamma(n, shape = 4, scale = 1 / 4)
  simDat <- reReg::simGSC(
    n = n, para = para, xmat = X, tau = 60,
    frailty = gamma, censoring = C, summary = TRUE
  )

  list(data = simDat, alpha_true = alpha.star, beta_true = beta.star)
}
