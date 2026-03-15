# Solve for Baseline Recurrence Hazard lambda0_r

Computes the closed-form solution for the baseline hazard point masses
of the recurrence sub-model.

## Usage

``` r
jfm_lambda0r_solution(tr, d_tr, n, S0t_re)
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

## Value

Named numeric vector of baseline hazard point masses for recurrence.
