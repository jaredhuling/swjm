# Compute the Cooperative Lasso Norm

For each pair (theta_j, theta\_(j+p)), computes the cooperative lasso
norm: L2 norm if signs agree, L1 norm if signs disagree.

## Usage

``` r
coop_norm(theta, p)
```

## Arguments

- theta:

  Numeric vector of length 2p (concatenation of alpha and beta).

- p:

  Integer, number of covariates.

## Value

Scalar cooperative lasso norm value.
