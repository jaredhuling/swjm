# Generate Both Baseline Hazards (JFM)

Computes Breslow-type baseline hazards for both death and recurrence.

## Usage

``` r
lambda0_gen(Data2, alpha, beta)
```

## Arguments

- Data2:

  A data frame in recurrent-event format.

- alpha:

  Coefficient vector for the recurrence sub-model (first p elements of
  theta).

- beta:

  Coefficient vector for the death sub-model (second p elements of
  theta).

## Value

A list with components `lambda0_r` and `lambda0_d`.
