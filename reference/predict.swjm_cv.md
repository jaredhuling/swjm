# Predict from a Cross-Validated Joint Model Fit

Computes subject-specific predictions from a cross-validated `swjm_cv`
fit. The output differs by model:

## Usage

``` r
# S3 method for class 'swjm_cv'
predict(object, newdata, times = NULL, ...)
```

## Arguments

- object:

  An object of class `"swjm_cv"`.

- newdata:

  A numeric matrix or data frame of covariate values. Must have `p`
  columns named `x1`, ..., `xp` (if a data frame) or exactly `p` columns
  (if a matrix).

- times:

  Numeric vector of evaluation times (JFM only). If `NULL`, the observed
  event times from the training data are used.

- ...:

  Currently unused.

## Value

An object of class `"swjm_pred"`, a list with:

- S_re:

  Matrix of readmission-free survival probabilities (rows = subjects,
  columns = `times`).

- S_de:

  Matrix of death survival probabilities.

- times:

  Numeric vector of evaluation times.

- lp_re:

  Linear predictors for readmission (\\\hat\alpha^\top z_i\\). For JSCM
  this is the log time-acceleration for the recurrence process.

- lp_de:

  Linear predictors for death (\\\hat\beta^\top z_i\\). For JSCM this is
  the log time-acceleration for the terminal process.

- time_accel_re:

  (JSCM only) \\e^{\hat\alpha^\top z_i}\\: the multiplicative factor by
  which the recurrence time axis is scaled relative to baseline. `NULL`
  for JFM.

- time_accel_de:

  (JSCM only) Analogous time-acceleration factor for the terminal
  process. `NULL` for JFM.

- contrib_re:

  Matrix of per-predictor contributions \\\hat\alpha_j z\_{ij}\\ (rows =
  subjects, columns = covariates). For JFM these are log-hazard
  contributions; for JSCM they are log time-acceleration contributions.

- contrib_de:

  Analogous matrix for the terminal process.

## Details

**JFM** — returns survival curves for both processes using the Breslow
cumulative baseline hazards. For subject \\i\\ with covariate vector
\\z_i\\: \$\$ S\_{\text{re}}(t \mid z_i) =
\exp\\\bigl(-\hat\Lambda_0^r(t)\\ e^{\hat\alpha^\top z_i}\bigr), \quad
S\_{\text{de}}(t \mid z_i) = \exp\\\bigl(-\hat\Lambda_0^d(t)\\
e^{\hat\beta^\top z_i}\bigr). \$\$

**JSCM** — returns survival curves for both processes using a
Nelson-Aalen baseline on the accelerated time scale. For subject \\i\\:
\$\$ S\_{\text{re}}(t \mid z_i) = \exp\\\bigl(-\hat\Lambda_0^r(t\\
e^{\hat\alpha^\top z_i})\bigr), \quad S\_{\text{de}}(t \mid z_i) =
\exp\\\bigl(-\hat\Lambda_0^d(t\\ e^{\hat\beta^\top z_i})\bigr). \$\$ The
linear predictor \\\hat\alpha^\top z_i\\ is a log time-acceleration
factor: the recurrence process for subject \\i\\ runs on a time axis
scaled by \\e^{\hat\alpha^\top z_i}\\ relative to the baseline. Each
predictor contributes \\\hat\alpha_j z\_{ij}\\ to this log-scale factor,
so \\e^{\hat\alpha_j z\_{ij}}\\ is the multiplicative factor on the time
scale attributable to predictor \\j\\ alone. Values greater than 1
accelerate events (shorter times); values less than 1 decelerate them
(longer times).

## Examples

``` r
# \donttest{
dat <- generate_data(n = 50, p = 5, scenario = 1, model = "jfm")
#> Error in beta %*% z: non-conformable arguments
cv  <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
                    max_iter = 100)
#> Error: object 'dat' not found
newz <- matrix(rnorm(15), nrow = 3, ncol = 5)
pred <- predict(cv, newdata = newz)
#> Error: object 'cv' not found
plot(pred)
#> Error: object 'pred' not found

dat_jscm <- generate_data(n = 50, p = 5, scenario = 1, model = "jscm")
#> Error in reReg::simGSC(n = n, para = para, xmat = X, tau = 60, frailty = gamma,     censoring = C, summary = TRUE): Parameter alpha does not match with the number of covariates.
cv_jscm  <- cv_stagewise(dat_jscm$data, model = "jscm", penalty = "coop",
                          max_iter = 500)
#> Error: object 'dat_jscm' not found
pred_jscm <- predict(cv_jscm, newdata = newz)
#> Error: object 'cv_jscm' not found
plot(pred_jscm)
#> Error: object 'pred_jscm' not found
# }
```
