# Solve for Baseline Death Hazard lambda0_d

Computes the closed-form solution for the baseline hazard point masses
of the death sub-model.

## Usage

``` r
jfm_lambda0d_solution(td, d_td, n, S0t_de)
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

## Value

Named numeric vector of baseline hazard point masses for death.
