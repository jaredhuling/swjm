# Estimating Equation for Beta (Death) in the JSCM

Evaluates the estimating equation for the terminal event coefficient
vector beta in the joint scale-change model, using the C++
implementations `temLog` and `reRate`.

## Usage

``` r
jscm_ee_beta(Data2, a, b)
```

## Arguments

- Data2:

  A data frame in recurrent-event format.

- a:

  Numeric vector of current alpha coefficients.

- b:

  Numeric vector of current beta coefficients.

## Value

Numeric vector of score values, one per covariate.
