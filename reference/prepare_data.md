# Prepare Data for swjm Functions

Expands factor and character covariate columns into dummy (indicator)
variables. The first 5 columns (id, t.start, t.stop, event, status) are
left unchanged. For each factor/character column with L levels, L-1
dummy columns are created (reference level dropped). A message is
printed when expansion occurs.

## Usage

``` r
prepare_data(data, caller = "swjm")
```

## Arguments

- data:

  A data frame.

- caller:

  Character string naming the calling function (for messages).

## Value

A data frame with factor/character covariates replaced by numeric dummy
columns.
