# Score Equation U2 for Alpha (Death)

Evaluates the score equation for the death coefficient vector alpha.

## Usage

``` r
jfm_score_alpha(pseudo_entries, index_death_matrix, td_id, S1t_de, S0t_de)
```

## Arguments

- pseudo_entries:

  Sorted pseudo data set.

- index_death_matrix:

  Matrix of pseudo entry indices for death times.

- td_id:

  Reordered subject IDs for death events.

- S1t_de:

  S1(t) matrix for death.

- S0t_de:

  S0(t) vector for death.

## Value

Numeric vector of score values, one per covariate.
