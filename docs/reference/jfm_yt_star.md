# At-Risk Indicator Y\*(t) with Status-Dependent Comparison

Returns 1 if the subject is at risk at time t, using strict or
non-strict inequality depending on the status indicator.

## Usage

``` r
jfm_yt_star(t, y, stat)
```

## Arguments

- t:

  A scalar time point.

- y:

  The subject's observed time.

- stat:

  The subject's status indicator (1 = event, 0 = censored).

## Value

Integer 0 or 1 indicator.
