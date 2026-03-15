# G-bar(t) Evaluated at Each Death Time

Computes the ratio G-bar(t) at each death event time, used in the
estimating equation for theta.

## Usage

``` r
jfm_gt_bar_death(r2i_death_subject_matrix, td, Y, STATUS, list_recur)
```

## Arguments

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

## Value

A matrix (1 row) of G-bar values, one per death time.
