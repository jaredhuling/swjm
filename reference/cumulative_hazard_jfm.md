# Cumulative Hazard from Estimated Baseline Hazards (JFM)

Computes subject-specific cumulative hazards for recurrent events using
estimated baseline hazard point masses and the pseudo data approach.

## Usage

``` r
cumulative_hazard_jfm(t.start, I, Z, beta, tr, lambda0_r, tr.id)
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

  Matrix of cumulative hazard values (rows = subjects, cols = recurrent
  times).

- index_recurrent_matrix:

  Matrix of pseudo entry indices.

- tr_id:

  Reordered subject IDs for recurrent events.
