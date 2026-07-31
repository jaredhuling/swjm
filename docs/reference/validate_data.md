# Validate Recurrent-Event Data Frame

Checks that a data frame has the required structure for swjm functions:
columns id, t.start, t.stop, event, status, and at least one covariate.
Produces informative error messages for each violated requirement.

## Usage

``` r
validate_data(data, caller = "swjm")
```

## Arguments

- data:

  Object to validate.

- caller:

  Character string naming the calling function (for error messages).

## Value

Invisibly returns `TRUE` if valid; throws an error otherwise.
