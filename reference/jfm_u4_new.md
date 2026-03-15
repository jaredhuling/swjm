# Estimating Equation U4 (New) for Numerically Solving Baseline Death Hazard

Self-contained version of U4 that rebuilds the weight matrix and S0(t)
internally, used for numerically solving for lambda0_d when a
closed-form solution is not applicable.

## Usage

``` r
jfm_u4_new(lambda0_d, d_td, n, Y, td, theta, alpha, t.start, I, Z, td.id)
```

## Arguments

- lambda0_d:

  Baseline hazard point masses for death.

- d_td:

  Table of death time frequencies.

- n:

  Number of subjects.

- Y:

  Vector of composite censoring/death times per subject.

- td:

  Vector of death event times.

- theta:

  Frailty variance parameter.

- alpha:

  Coefficient vector for the death sub-model.

- t.start:

  Vector of interval start times.

- I:

  Vector of subject indicators for each pseudo entry.

- Z:

  List of covariate matrices, one per subject.

- td.id:

  Subject IDs corresponding to each death event.

## Value

Named numeric vector of U4 values, one per death time.
