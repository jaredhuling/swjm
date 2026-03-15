# Generate a Single Covariate Vector

Generates n observations from a specified distribution for simulation
purposes.

## Usage

``` r
jfm_v_gene(n, v.type, v.coeffi)
```

## Arguments

- n:

  Number of observations.

- v.type:

  Distribution type: 1 = normal, 2 = Bernoulli, 3 = Poisson.

- v.coeffi:

  Numeric vector of distribution parameters. For normal: c(mean, sd).
  For Bernoulli: c(prob). For Poisson: c(lambda).

## Value

Numeric vector of length n.
