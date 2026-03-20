# Cause-Specific CIF Markers (JFM)

Computes cause-specific cumulative incidence function (CIF) markers for
recurrent events under the joint frailty model, using estimated baseline
hazards.

## Usage

``` r
markers_cause_specific(alpha, beta, tr, td, lambda0_r, lambda0_d, Z_base)
```

## Arguments

- alpha:

  Coefficient vector for the recurrence sub-model (first p elements of
  theta).

- beta:

  Coefficient vector for the death sub-model (second p elements of
  theta).

- tr:

  Vector of recurrent event times.

- td:

  Vector of death event times.

- lambda0_r:

  Baseline hazard point masses for recurrence.

- lambda0_d:

  Baseline hazard point masses for death.

- Z_base:

  Matrix of baseline covariate values (one row per subject).

## Value

A list with components:

- cif_marker:

  Matrix of CIF values (rows = subjects, cols = event times).

- linear_marker:

  Numeric vector of linear predictors for recurrence.
