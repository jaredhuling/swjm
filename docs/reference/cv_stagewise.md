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
  pp = NULL
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

  Integer. Maximum iterations (passed to `stagewise_fit`).

- pp:

  Integer. Early-stop block size (passed to `stagewise_fit`).

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
#>   Best position (combined):    101  (lambda = 1.051)
#>   Selected variables:          3 readmission, 3 death
#>     Readmission (alpha): 1, 9, 10
#>     Death (beta):        1, 9, 10
# }
```
