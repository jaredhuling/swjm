# S1(t) for Recurrent Event Sub-Model with Cross-Fitting

Computes the weighted first moment S1(t) at each recurrent event time
for the recurrence sub-model, using fold-specific beta coefficients for
cross-fitting.

## Usage

``` r
jfm_s1t_recurrent_cf(
  Y,
  wt_recurrent_subject,
  tr,
  index_recurrent_matrix,
  pseudo_entries,
  beta_mat,
  CV_map
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

- beta_mat:

  Matrix of beta coefficients (rows = folds).

- CV_map:

  Two-column matrix mapping fold index to subject index.

## Value

Matrix of S1(t) values (rows = covariates, cols = recurrent times).
