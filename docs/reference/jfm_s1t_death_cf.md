# S1(t) for Death Sub-Model with Cross-Fitting

Computes the weighted first moment S1(t) at each death event time for
the death sub-model, using fold-specific alpha coefficients for
cross-fitting.

## Usage

``` r
jfm_s1t_death_cf(
  Y,
  wt_matrix,
  td,
  index_death_matrix,
  pseudo_entries,
  alpha_mat,
  CV_map
)
```

## Arguments

- Y:

  Vector of composite censoring/death times per subject.

- wt_matrix:

  Weight matrix from `jfm_wt_death`.

- td:

  Vector of death event times.

- index_death_matrix:

  Matrix of pseudo entry indices for death times.

- pseudo_entries:

  Sorted pseudo data set.

- alpha_mat:

  Matrix of alpha coefficients (rows = folds).

- CV_map:

  Two-column matrix mapping fold index to subject index.

## Value

Matrix of S1(t) values (rows = covariates, cols = death times).
