# S0(t) for Death Sub-Model

Computes the weighted zeroth moment S0(t) at each death event time for
the death sub-model.

## Usage

``` r
jfm_s0t_death(Y, wt_matrix, td, index_death_matrix, pseudo_entries, alpha)
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

- alpha:

  Coefficient vector for the death sub-model.

## Value

Numeric vector of S0(t) values, one per death time.
