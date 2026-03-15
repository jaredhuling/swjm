# Interpolate Coefficient Vectors Along a Lambda Sequence

Given a decreasing lambda sequence from a fitted path and a new set of
lambda values, linearly interpolates the corresponding coefficient
vectors.

## Usage

``` r
lambda_interp(lambda, s)
```

## Arguments

- lambda:

  Numeric vector of decreasing lambda values from the fitted path.

- s:

  Numeric vector of new lambda values at which to interpolate.

## Value

A list with components:

- left:

  Integer vector of left neighbor indices.

- right:

  Integer vector of right neighbor indices.

- frac:

  Numeric vector of interpolation fractions. The interpolated value is
  `frac * path[left] + (1 - frac) * path[right]`.
