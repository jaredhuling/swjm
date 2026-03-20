# Theta Estimating Equation for the Stagewise Algorithm

Self-contained version of the theta estimating equation that rebuilds
all intermediate quantities (wt_matrix, r2i, G-bar) internally, used as
the objective in a root-finding procedure for theta.

## Usage

``` r
jfm_est_theta_new(
  theta,
  t.start,
  I,
  Z,
  alpha,
  beta,
  lambda0_r,
  lambda0_d,
  td,
  tr,
  tr.id,
  td.id,
  Y,
  STATUS,
  list_recur,
  num_recur
)
```

## Arguments

- theta:

  Frailty variance parameter.

- t.start:

  Vector of interval start times.

- I:

  Vector of subject indicators for each pseudo entry.

- Z:

  List of covariate matrices, one per subject.

- alpha:

  Coefficient vector for the death sub-model.

- beta:

  Coefficient vector for the recurrent event sub-model.

- lambda0_r:

  Baseline hazard point masses for recurrence.

- lambda0_d:

  Baseline hazard point masses for death.

- td:

  Vector of death event times.

- tr:

  Vector of recurrent event times.

- tr.id:

  Subject IDs corresponding to each recurrent event.

- td.id:

  Subject IDs corresponding to each death event.

- Y:

  Vector of composite censoring/death times per subject.

- STATUS:

  Vector of death indicators per subject.

- list_recur:

  List of recurrent event times per subject.

- num_recur:

  Integer vector of recurrent event counts per subject.

## Value

Scalar value of the scaled estimating equation for theta.
