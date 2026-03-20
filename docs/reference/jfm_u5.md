# Estimating Equation U5 for Baseline Recurrence Hazard

Evaluates the estimating equation for the baseline hazard of the
recurrence sub-model at each recurrent event time.

## Usage

``` r
jfm_u5(tr, d_tr, n, S0t_re, lambda0_r)
```

## Arguments

- tr:

  Vector of recurrent event times.

- d_tr:

  Table of recurrent event frequencies.

- n:

  Number of subjects.

- S0t_re:

  S0(t) vector for recurrence.

- lambda0_r:

  Baseline hazard point masses for recurrence.

## Value

Named numeric vector of U5 values, one per recurrent event time.
