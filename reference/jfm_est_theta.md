# Estimating Equation for Theta

Evaluates the estimating equation for the frailty variance parameter
theta.

## Usage

``` r
jfm_est_theta(
  td_id,
  r2i_death_subject_matrix,
  td,
  Y,
  STATUS,
  list_recur,
  theta,
  num_recur
)
```

## Arguments

- td_id:

  Reordered subject IDs for death events.

- r2i_death_subject_matrix:

  Matrix of R2i values at death times.

- td:

  Vector of death event times.

- Y:

  Vector of composite censoring/death times per subject.

- STATUS:

  Vector of death indicators per subject.

- list_recur:

  List of recurrent event times per subject.

- theta:

  Current frailty variance parameter.

- num_recur:

  Integer vector of recurrent event counts per subject.

## Value

Scalar value of the estimating equation for theta.
