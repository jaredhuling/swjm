# Cause-Specific CIF Markers Under True Parameters (JFM)

Computes the CIF at a single time point under constant baseline hazards.

## Usage

``` r
markers_cause_specific_true(t, alpha, beta, lambda0_r, lambda0_d, Z_base)
```

## Arguments

- t:

  Scalar time at which to evaluate.

- alpha:

  Coefficient vector for the recurrence sub-model (first p elements of
  theta).

- beta:

  Coefficient vector for the death sub-model (second p elements of
  theta).

- lambda0_r:

  Scalar baseline hazard rate for recurrence.

- lambda0_d:

  Scalar baseline hazard rate for death.

- Z_base:

  Matrix of baseline covariate values (one row per subject).

## Value

A list with components:

- cif_marker:

  Numeric vector of CIF values per subject.

- linear_marker:

  Numeric vector of linear predictors for recurrence.
