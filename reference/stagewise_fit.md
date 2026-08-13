# Fit a Stagewise Regularization Path

Unified interface for stagewise variable selection in joint models of
recurrent events and terminal events. Dispatches to model-specific
implementations for the Joint Frailty Model (JFM) or Joint Scale-Change
Model (JSCM).

## Usage

``` r
stagewise_fit(
  data,
  model = c("jfm", "jscm"),
  penalty = c("coop", "lasso", "group"),
  eps = NULL,
  max_iter = NULL,
  pp = NULL,
  estimate_frailty = FALSE,
  standardize = TRUE,
  lambda_min_ratio = NULL
)
```

## Arguments

- data:

  A data frame in recurrent-event format with columns `id`, `t.start`,
  `t.stop`, `event`, `status`, and covariate columns `x1`, ..., `xp`.

- model:

  Character. Either `"jfm"` or `"jscm"`.

- penalty:

  Character. One of `"coop"` (cooperative lasso), `"lasso"`, or
  `"group"` (group lasso).

- eps:

  Numeric. Base step size for the stagewise update, used at the top of
  the regularization path. The step is divided by ten each time the dual
  norm falls a further decade, so `eps` sets the resolution of the whole
  path. If `NULL`, defaults to 0.1 for the JFM and 0.01 for the JSCM.
  Smaller values trace the path more finely at proportionally greater
  cost; `max_iter` must be large enough for the path to reach small
  `lambda`.

- max_iter:

  Integer. Maximum number of stagewise iterations.

- pp:

  Integer. Early-stopping block size: algorithm checks every `pp`
  iterations if fewer than 3 unique coordinates were updated.

- estimate_frailty:

  Logical. For JFM only: if `TRUE`, estimates the frailty variance and
  uses the Kalbfleisch et al. (2013) frailty weights \\w_i(t)\\ in the
  estimating equations. If `FALSE` (default), uses unit weights
  (simplified model without frailty).

- standardize:

  Logical. If `TRUE` (default), covariates are divided by their standard
  deviations before fitting and coefficients are rescaled back to the
  original scale. For the JFM with time-varying covariates, the SD is
  computed across all rows (all time points and subjects). For the JSCM
  with time-invariant covariates, the SD is computed across subjects
  only.

- lambda_min_ratio:

  Numeric. The path stops once the dual norm falls below this fraction
  of its value at the top of the path. Past that point the coefficients
  are essentially static, so the remaining iterations cost time without
  changing the fit. If `NULL`, defaults to 0.01 when `n > p` and 1e-4
  otherwise, following the convention used by `glmnet`. Set to 0 to
  trace the path to the smallest reachable `lambda`.

## Value

An object of class `"swjm_path"`, a list with components:

- k:

  Number of iterations performed.

- stop_reason:

  Character, why the path terminated: `"max_iter"` (ran to the iteration
  cap), `"early_stop"` (a single coordinate dominated every update in
  the last `pp` iterations), or `"lambda_min_ratio"` (the dual norm fell
  below `lambda_min_ratio` times its initial value).

- max_iter:

  Integer, the iteration cap used for the fit.

- theta:

  Matrix of coefficient paths (`2p` rows by `k+1` columns).

- lambda:

  Numeric vector of penalty parameter approximations at each iteration.

- norm:

  Numeric vector of penalty norm values along the path.

- model:

  Character, the model used.

- penalty:

  Character, the penalty used.

## Examples

``` r
# \donttest{
dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
                     max_iter = 100)
fit
#> Stagewise path (jfm/coop)
#> 
#>   Covariates (p):            10
#>   Iterations:                100
#>   Lambda range:              [0.6026, 1.363]
#>   Active at final step:      4 readmission, 4 death
#>     Readmission (alpha): 1, 2, 9, 10
#>     Death (beta):        1, 3, 9, 10
# }
```
