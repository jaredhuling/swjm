# Estimating Equation U4 for Baseline Death Hazard

Evaluates the estimating equation for the baseline hazard of the death
sub-model at each death event time.

## Usage

``` r
jfm_u4(td, d_td, n, S0t_de, lambda0_d)
```

## Arguments

- td:

  Vector of death event times.

- d_td:

  Table of death time frequencies.

- n:

  Number of subjects.

- S0t_de:

  S0(t) vector for death.

- lambda0_d:

  Baseline hazard point masses for death.

## Value

Named numeric vector of U4 values, one per death time.
