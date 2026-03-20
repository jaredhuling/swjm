# W_i(t) Evaluated at Each Death Time for Each Subject

Constructs the weight matrix W_i(t) evaluated at each death event time
for each subject, using the pseudo data set approach.

## Usage

``` r
jfm_wt_death(theta, alpha, t.start, I, Z, td, lambda0_d, td.id)
```

## Arguments

- theta:

  Frailty variance parameter.

- alpha:

  Coefficient vector for the death sub-model.

- t.start:

  Vector of interval start times.

- I:

  Vector of subject indicators for each pseudo entry.

- Z:

  List of covariate matrices, one per subject.

- td:

  Vector of death event times.

- lambda0_d:

  Baseline hazard point masses for death.

- td.id:

  Subject IDs corresponding to each death event.

## Value

A list with components:

- wt_matrix:

  Matrix of W_i(t) values (rows = subjects, cols = death times).

- td_id:

  Reordered subject IDs for death events.

- index_death_matrix:

  Matrix of pseudo entry indices (rows = subjects, cols = death times).

- pseudo_entries:

  Sorted pseudo data set.
