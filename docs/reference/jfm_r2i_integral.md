# Integral in R2i(t) Evaluated at Each Recurrent Event Time

Computes the integral component of R2i(t) for each subject at each
recurrent event time, using the pseudo data set approach.

## Usage

``` r
jfm_r2i_integral(t.start, I, Z, beta, tr, lambda0_r, tr.id)
```

## Arguments

- t.start:

  Vector of interval start times.

- I:

  Vector of subject indicators for each pseudo entry.

- Z:

  List of covariate matrices, one per subject.

- beta:

  Coefficient vector for the recurrent event sub-model.

- tr:

  Vector of recurrent event times.

- lambda0_r:

  Baseline hazard point masses for recurrence.

- tr.id:

  Subject IDs corresponding to each recurrent event.

## Value

A list with components:

- integral_matrix:

  Matrix of integral values (rows = subjects, cols = recurrent times).

- index_recurrent_matrix:

  Matrix of pseudo entry indices.

- tr_id:

  Reordered subject IDs for recurrent events.
