# Extract Decreasing Lambda Path

Post-processes a raw lambda sequence from the stagewise algorithm to
keep only strictly decreasing values (via running minimum and
deduplication).

## Usage

``` r
extract_decreasing_indices(lambda, tol_digits = 6L)
```

## Arguments

- lambda:

  Numeric vector of raw lambda values from stagewise iterations.

- tol_digits:

  Integer, number of digits for rounding in deduplication.

## Value

An integer vector of indices into the original lambda vector
corresponding to the strictly decreasing subsequence.
