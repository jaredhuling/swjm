# Score Equation U1 for Beta (Recurrence)

Evaluates the score equation for the recurrence coefficient vector beta.

## Usage

``` r
jfm_score_beta(pseudo_entries, index_recurrent_matrix, tr_id, S1t_re, S0t_re)
```

## Arguments

- pseudo_entries:

  Sorted pseudo data set.

- index_recurrent_matrix:

  Matrix of pseudo entry indices for recurrent times.

- tr_id:

  Reordered subject IDs for recurrent events.

- S1t_re:

  S1(t) matrix for recurrence.

- S0t_re:

  S0(t) vector for recurrence.

## Value

Numeric vector of score values, one per covariate.
