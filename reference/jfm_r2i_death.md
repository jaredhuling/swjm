# R2i(t) Evaluated at Each Death Time for Each Subject

Computes R2i(t) for each subject at each death event time, by combining
the recurrent event integral with the weight matrix.

## Usage

``` r
jfm_r2i_death(t.start, I, Z, beta, tr, lambda0_r, td, wt_matrix, tr.id)
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

- td:

  Vector of death event times.

- wt_matrix:

  Weight matrix from `jfm_wt_death`.

- tr.id:

  Subject IDs corresponding to each recurrent event.

## Value

A list with components:

- r2i_death_subject_matrix:

  Matrix of R2i values at death times.

- index_recurrent_matrix:

  Matrix of pseudo entry indices.

- tr_id:

  Reordered subject IDs for recurrent events.
