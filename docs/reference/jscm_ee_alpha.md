# Estimating Equation for Alpha (Recurrence) in the JSCM

Evaluates the estimating equation for the recurrence coefficient vector
alpha in the joint scale-change model, using the C++ implementation
`am1`.

## Usage

``` r
jscm_ee_alpha(Data2, a)
```

## Arguments

- Data2:

  A data frame in recurrent-event format.

- a:

  Numeric vector of current alpha coefficients.

## Value

Numeric vector of score values, one per covariate.
