# Extract Common Data Components from a Recurrent-Event Data Frame

Parses a standard recurrent-event data frame into the components needed
by the JFM estimating equations: covariate lists, event times, at-risk
indicators, etc.

## Usage

``` r
extract_data_components(Data2)
```

## Arguments

- Data2:

  A data frame in the standard recurrent-event format.

## Value

A list with components:

- Z:

  List of covariate matrices, one per subject.

- n:

  Number of unique subjects.

- p:

  Number of covariates.

- td:

  Death times.

- td.id:

  Subject IDs for death times.

- d_td:

  Table of death time frequencies.

- tr:

  Recurrent event times.

- tr.id:

  Subject IDs for recurrent events.

- d_tr:

  Table of recurrent event frequencies.

- Y:

  Composite censoring/death times (one per subject).

- STATUS:

  Death indicator at composite time (one per subject).

- list_recur:

  List of recurrent event times per subject.

- num_recur:

  Integer vector of recurrent event counts per subject.

- t.start:

  All interval start times.

- I:

  All subject IDs (matching rows of the data frame).
