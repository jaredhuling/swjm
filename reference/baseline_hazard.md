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
dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
cv  <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
                    max_iter = 100)
bh  <- baseline_hazard(cv)
head(bh)
#>          time cumhaz_readmission cumhaz_death
#> 1 0.002991397         0.01161349            0
#> 2 0.005471293         0.02348072            0
#> 3 0.007580442         0.03593013            0
#> 4 0.010431138         0.04950762            0
#> 5 0.015063521         0.06431059            0
#> 6 0.016583838         0.07976426            0
# }
```
