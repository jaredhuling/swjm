# Cross-Validation for Stagewise Variable Selection

Selects the optimal penalty parameter (lambda) along the stagewise path
using K-fold cross-validation with cross-fitted estimating equations.

## Usage

``` r
cv_stagewise(
  data,
  model = c("jfm", "jscm"),
  penalty = c("coop", "lasso", "group"),
  K = 5L,
  lambda_seq = NULL,
  eps = NULL,
  max_iter = NULL,
  pp = NULL,
  estimate_frailty = FALSE,
  ncores = 1L,
  standardize = TRUE,
  lambda_min_ratio = NULL
)
```

## Arguments

- data:

  A data frame in recurrent-event format.

- model:

  Character. Either `"jfm"` or `"jscm"`.

- penalty:

  Character. One of `"coop"`, `"lasso"`, or `"group"`.

- K:

  Integer. Number of cross-validation folds (default 5).

- lambda_seq:

  Numeric vector. The lambda sequence at which to evaluate the
  cross-validation criterion. If `NULL`, extracted from a full-data
  stagewise fit.

- eps:

  Numeric. Step size (passed to `stagewise_fit`).

- max_iter:

  Integer. Maximum iterations (passed to `stagewise_fit`). A fold path
  cut off by this cap before reaching the smallest lambda in
  `lambda_seq` is automatically refit with a doubled budget (up to 8
  times `max_iter`), so out-of-fold scores at the tail of the grid come
  from interpolation on a completed path rather than extrapolation from
  a truncated one; a warning is issued if coverage still fails.

- pp:

  Integer. Early-stop block size (passed to `stagewise_fit`).

- estimate_frailty:

  Logical. For JFM only: if `TRUE`, estimates the frailty variance and
  uses frailty weights in the estimating equations (passed to
  `stagewise_fit`).

- ncores:

  Integer. Number of cores for parallel fold training (default 1,
  sequential). Uses
  [`parallel::parLapply`](https://rdrr.io/r/parallel/clusterApply.html)
  with a PSOCK cluster, which works on all platforms including Windows.

- standardize:

  Logical. If `TRUE` (default), covariates are standardized before
  fitting (passed to `stagewise_fit`).

- lambda_min_ratio:

  Numeric. Passed to
  [`stagewise_fit`](http://jaredhuling.org/swjm/reference/stagewise_fit.md)
  for the full-data path; see there for the default. Fold fits are
  stopped at the smallest lambda in `lambda_seq` regardless (with
  iteration-budget doubling as described under `max_iter`), so they
  cover the grid being scored.

## Value

An object of class `"swjm_cv"`, a list with components:

- position_CF:

  Integer, position of best lambda by combined cross-fitted score norm.

- position_CF_re:

  Integer, position of best lambda by recurrence score norm.

- position_CF_cen:

  Integer, position of best lambda by terminal score norm.

- lambda_seq:

  Numeric vector of lambda values evaluated.

- Scorenorm_crossfit:

  Combined cross-fitted score norm path.

- Scorenorm_crossfit_re:

  Recurrence score norm path.

- Scorenorm_crossfit_ce:

  Terminal score norm path.

- full_fit:

  The full-data stagewise fit (class `"swjm_path"`).

- model:

  Character.

- penalty:

  Character.

## Examples

``` r
# \donttest{
dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
                     max_iter = 100)
cv_res <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
                       lambda_seq = fit$lambda, max_iter = 100)
cv_res
#> Cross-validation (jfm/coop)
#> 
#>   Covariates (p):              10
#>   Lambda grid size:            101
#>   Best position (combined):    101  (lambda = 0.6182)
#>   Selected variables:          5 readmission, 6 death
#>     Readmission (alpha): 1, 2, 4, 9, 10
#>     Death (beta):        1, 2, 3, 4, 9, 10
# }
```
