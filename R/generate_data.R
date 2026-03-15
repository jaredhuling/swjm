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
#'   For JFM: \code{b}, \code{lambda0_d}, \code{lambda0_r}.
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
                               lambda0_d = 0.041, lambda0_r = 1) {
  S <- scenario
  theta <- 1
  alpha.star <- numeric(p)
  beta.star  <- numeric(p)

  if (S == 1) {
    # first 10% large alpha and small beta (+)
    beta.star[1:(10 * 0.1)] <- 1.1
    alpha.star[1:(10 * 0.1)] <- 0.1
    # 10%-20% large alpha and small beta (-)
    beta.star[(10 * 0.1 + 1):(10 * 0.2)] <- -1.1
    alpha.star[(10 * 0.1 + 1):(10 * 0.2)] <- -0.1
    # 20-30% small alpha and large beta (+)
    beta.star[(10 * 0.2 + 1):(10 * 0.3)] <- 0.1
    alpha.star[(10 * 0.2 + 1):(10 * 0.3)] <- 1.1
    # 30%-40% small alpha and large beta (-)
    beta.star[(10 * 0.3 + 1):(10 * 0.4)] <- -0.1
    alpha.star[(10 * 0.3 + 1):(10 * 0.4)] <- -1.1
    # 80-90% similar magnitude (+)
    beta.star[(10 * 0.8 + 1):(10 * 0.9)] <- 1
    alpha.star[(10 * 0.8 + 1):(10 * 0.9)] <- 1
    # 90-100% similar magnitude (-)
    beta.star[(10 * 0.9 + 1):(10 * 1)] <- -1
    alpha.star[(10 * 0.9 + 1):(10 * 1)] <- -1

  } else if (S == 2) {
    # give the 10% to 20% both strong signal and 20% to 30% both weak signal
    beta.star[(10 * 0.1 + 1):(10 * 0.2)] <- -1
    beta.star[(10 * 0.2 + 1):(10 * 0.3)] <- -1

    # first 10%
    beta.star[1:(10 * 0.1)] <- 1

    # 30% to 40%
    beta.star[(10 * 0.3 + 1):(10 * 0.4)] <- 1

    # 80% to 90%
    beta.star[(10 * 0.8 + 1):(10 * 0.9)] <- 1

    # 0-10%
    alpha.star[1:(10 * 0.1)] <- 1
    # 10%-20%
    alpha.star[(10 * 0.1 + 1):(10 * 0.2)] <- -1
    # 20%-30%
    alpha.star[(10 * 0.2 + 1):(10 * 0.3)] <- 1
    # 30%-40%
    alpha.star[(10 * 0.3 + 1):(10 * 0.4)] <- -1
    # 90% to 100%
    alpha.star[(10 * 0.9 + 1):10] <- 1

  } else if (S == 3) {
    # 10%: alpha:1 beta: -1
    beta.star[1:(10 * 0.1)] <- 1
    alpha.star[1:(10 * 0.1)] <- -1
    # 10-20%: alpha:-1 beta: 1
    beta.star[(10 * 0.1 + 1):(10 * 0.2)] <- -1
    alpha.star[(10 * 0.1 + 1):(10 * 0.2)] <- 1
    # 20-30%: alpha:2 beta: -1
    beta.star[(10 * 0.2 + 1):(10 * 0.3)] <- 1.5
    alpha.star[(10 * 0.2 + 1):(10 * 0.3)] <- -1
    # 30-40%: alpha:-1 beta: 2
    beta.star[(10 * 0.3 + 1):(10 * 0.4)] <- -1
    alpha.star[(10 * 0.3 + 1):(10 * 0.4)] <- 1.5
    # 80-90%: alpha:1 beta: 0
    beta.star[(10 * 0.8 + 1):(10 * 0.9)] <- 1
    alpha.star[(10 * 0.8 + 1):(10 * 0.9)] <- 0
    # 90-100%: alpha:0 beta: 1
    beta.star[(10 * 0.9 + 1):10] <- 0
    alpha.star[(10 * 0.9 + 1):10] <- 1
  } else {
    beta.star <- alpha.star <- c(-0.5, 0.5)
  }

  alpha <- alpha.star
  beta  <- beta.star

  id         <- c()
  t.start    <- c()
  t.stop     <- c()
  event      <- c()
  status     <- c()
  cov.matrix <- c()

  for (j in 1:n) {
    censor <- runif(1, 1, b)
    lam_d  <- lambda0_d
    lam_r  <- lambda0_r
    gamma  <- 1
    Xk     <- 0
    delta  <- 0

    X.vec      <- c()
    sub.event  <- c()
    sub.status <- c()
    cov.vec    <- c()
    sub.id     <- c()

    while (Xk < censor & delta == 0) {
      z       <- rnorm(p)
      cov.vec <- rbind(cov.vec, z)

      T_gap <- rexp(1, rate = gamma * lam_r * exp(drop(beta %*% z)))
      D_gap <- rexp(1, rate = gamma * lam_d * exp(drop(alpha %*% z)))

      Xk    <- Xk + min(T_gap, D_gap)
      X.vec <- c(X.vec, Xk)

      if (D_gap < T_gap & Xk < censor) {
        delta <- 1
      } else {
        delta <- 0
      }

      sub.status <- c(sub.status, delta)
      sub.event  <- c(sub.event, 1 - delta)
      sub.id     <- c(sub.id, j)
    }

    if (delta != 1) {
      X.vec[length(X.vec)]         <- censor
      sub.event[length(sub.event)] <- 0
    }

    sub.t.stop  <- X.vec
    sub.t.start <- c(0, X.vec[-length(X.vec)])

    id         <- c(id, sub.id)
    t.start    <- c(t.start, sub.t.start)
    t.stop     <- c(t.stop, sub.t.stop)
    event      <- c(event, sub.event)
    status     <- c(status, sub.status)
    cov.matrix <- rbind(cov.matrix, cov.vec)
  }

  colnames(cov.matrix) <- paste0("x", 1:p)
  Data2 <- data.frame(
    id = id, t.start = t.start, t.stop = t.stop,
    event = event, status = status, cov.matrix
  )
  rownames(Data2) <- 1:nrow(Data2)

  # Relabel to package convention: alpha = recurrence, beta = death
  # (internally alpha.star drove death and beta.star drove recurrence)
  list(data = Data2, alpha_true = beta.star, beta_true = alpha.star)
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
