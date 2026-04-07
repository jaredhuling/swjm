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
  standardize = TRUE
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

  Numeric. Step size for the stagewise update. If `NULL`, uses adaptive
  step size.

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

## Value

An object of class `"swjm_path"`, a list with components:

- k:

  Number of iterations performed.

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
#>   Lambda range:              [1.069, 1.363]
#>   Active at final step:      4 readmission, 3 death
#>     Readmission (alpha): 1, 2, 9, 10
#>     Death (beta):        1, 9, 10
# }
```
