# Cumulative Hazard Under True Baseline Hazard (JFM)

Computes the cumulative hazard at time `x` under piecewise-constant
baseline hazard with known cut points and covariate values.

## Usage

``` r
cumulative_hazard_true_jfm(x, cut_points, z_values, beta)
```

## Arguments

- x:

  Scalar time at which to evaluate.

- cut_points:

  Numeric vector of cut point times.

- z_values:

  Matrix of covariate values per interval.

- beta:

  Coefficient vector.

## Value

Scalar cumulative hazard value.
