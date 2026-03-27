# Generate Simulated Data for the Joint Scale-Change Model (JSCM)

Generates recurrent-event and terminal-event data under an AFT-type
joint scale-change model using
[`simGSC`](https://rdrr.io/pkg/reReg/man/simGSC.html). Ported from
`Data_gen_reReg()`.

## Usage

``` r
generate_data_jscm(n, p, scenario = 1, b = 4)
```

## Arguments

- n:

  Integer. Number of subjects.

- p:

  Integer. Number of covariates.

- scenario:

  Integer. Scenario (1, 2, 3, or other).

- b:

  Numeric. Upper bound of the censoring uniform distribution (default
  4).

## Value

A list with components:

- data:

  Object returned by
  [`simGSC`](https://rdrr.io/pkg/reReg/man/simGSC.html) (a data frame
  with recurrent-event structure).

- alpha_true:

  True alpha (recurrence) coefficients.

- beta_true:

  True beta (terminal) coefficients.

## Details

In the JSCM convention:

- `alpha` governs the recurrence process.

- `beta` governs the terminal (death) process.

Covariates are drawn from `Uniform(-1, 1)`. A gamma frailty with shape =
4, scale = 1/4 is used. Censoring times are `Uniform(0, b)`.

## Examples

``` r
# \donttest{
dat <- generate_data_jscm(n = 30, p = 10, scenario = 1)
#> Call: 
#> reReg::simGSC(n = n, summary = TRUE, para = para, xmat = X, censoring = C, 
#>     frailty = gamma, tau = 60)
#> 
#> Summary:
#> Sample size:                                    30 
#> Number of recurrent event observed:             65 
#> Average number of recurrent event per subject:  2.167 
#> Proportion of subjects with a terminal event:   0.167 
#> 
#> 
head(dat$data)
#>   id    t.start     t.stop event status         x1        x2        x3
#> 1  1 0.00000000 0.09560175     1      0 -0.3892932 0.3198084 0.4810566
#> 2  1 0.09560175 0.16089864     1      0 -0.3892932 0.3198084 0.4810566
#> 3  1 0.16089864 0.64952068     1      0 -0.3892932 0.3198084 0.4810566
#> 4  1 0.64952068 0.70419923     1      0 -0.3892932 0.3198084 0.4810566
#> 5  1 0.70419923 0.79181158     1      0 -0.3892932 0.3198084 0.4810566
#> 6  1 0.79181158 0.87123659     1      0 -0.3892932 0.3198084 0.4810566
#>           x4          x5       x6         x7        x8        x9        x10
#> 1 -0.4277874 -0.02930706 0.746448 -0.1607343 0.6105984 0.9964618 -0.7943127
#> 2 -0.4277874 -0.02930706 0.746448 -0.1607343 0.6105984 0.9964618 -0.7943127
#> 3 -0.4277874 -0.02930706 0.746448 -0.1607343 0.6105984 0.9964618 -0.7943127
#> 4 -0.4277874 -0.02930706 0.746448 -0.1607343 0.6105984 0.9964618 -0.7943127
#> 5 -0.4277874 -0.02930706 0.746448 -0.1607343 0.6105984 0.9964618 -0.7943127
#> 6 -0.4277874 -0.02930706 0.746448 -0.1607343 0.6105984 0.9964618 -0.7943127
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
# }
```
