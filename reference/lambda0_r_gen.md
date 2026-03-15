# Generate Baseline Hazard for Recurrence (JFM)

Computes the Breslow-type baseline hazard for recurrence from a fitted
alpha coefficient vector.

## Usage

``` r
lambda0_r_gen(Data2, alpha)
```

## Arguments

- Data2:

  A data frame in recurrent-event format.

- alpha:

  Coefficient vector for the recurrence sub-model (first p elements of
  theta).

## Value

Named numeric vector of baseline hazard point masses.
