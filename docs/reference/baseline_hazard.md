# Cumulative Baseline Hazard Step Functions

Evaluates the cumulative baseline hazards for readmission and death at
arbitrary time points. For JFM, Breslow-type estimators are used. For
JSCM, Nelson-Aalen estimators on the accelerated time scale are used.

## Usage

``` r
baseline_hazard(
  object,
  times = NULL,
  which = c("both", "readmission", "death")
)
```

## Arguments

- object:

  An object of class `"swjm_cv"`.

- times:

  Numeric vector of evaluation times. If `NULL`, the observed event
  times stored in the fit are used.

- which:

  Character. Which baseline hazard(s) to return: `"both"` (default),
  `"readmission"`, or `"death"`.

## Value

A data frame with column `time` and, depending on `which`,
`cumhaz_readmission` and/or `cumhaz_death`.

## Examples

``` r
# \donttest{
dat <- generate_data(n = 50, p = 5, scenario = 1, model = "jfm")
#> Error in beta %*% z: non-conformable arguments
cv  <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
                    max_iter = 100)
#> Error: object 'dat' not found
bh  <- baseline_hazard(cv)
#> Error: object 'cv' not found
head(bh)
#> Error: object 'bh' not found
# }
```
