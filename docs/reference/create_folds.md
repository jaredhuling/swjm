# Create Stratified K-Fold Splits

Randomly assigns a vector of IDs into K approximately equal-sized folds.

## Usage

``` r
create_folds(ids, K)
```

## Arguments

- ids:

  Vector of subject IDs.

- K:

  Integer number of folds.

## Value

A list of length K, where each element contains the IDs assigned to that
fold.
