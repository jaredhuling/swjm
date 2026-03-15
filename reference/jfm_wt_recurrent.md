# W_i(t) Evaluated at Each Recurrent Event Time for Each Subject

Constructs the weight matrix W_i(t) evaluated at each recurrent event
time for each subject, by mapping recurrent event times to the nearest
death time in the weight matrix.

## Usage

``` r
jfm_wt_recurrent(tr, wt_matrix, td)
```

## Arguments

- tr:

  Vector of recurrent event times.

- wt_matrix:

  Weight matrix from `jfm_wt_death`.

- td:

  Vector of death event times.

## Value

Matrix of W_i(t) values at recurrent event times (rows = subjects, cols
= recurrent event times).
