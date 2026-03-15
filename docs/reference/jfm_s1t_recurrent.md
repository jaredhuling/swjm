# S1(t) for Recurrent Event Sub-Model

Computes the weighted first moment S1(t) at each recurrent event time
for the recurrence sub-model.

## Usage

``` r
jfm_s1t_recurrent(
  Y,
  wt_recurrent_subject,
  tr,
  index_recurrent_matrix,
  pseudo_entries,
  beta
)
```

## Arguments

- Y:

  Vector of composite censoring/death times per subject.

- wt_recurrent_subject:

  Weight matrix at recurrent event times.

- tr:

  Vector of recurrent event times.

- index_recurrent_matrix:

  Matrix of pseudo entry indices for recurrent times.

- pseudo_entries:

  Sorted pseudo data set.

- beta:

  Coefficient vector for the recurrent event sub-model.

## Value

Matrix of S1(t) values (rows = covariates, cols = recurrent times).
