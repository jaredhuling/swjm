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
#> Number of recurrent event observed:             63 
#> Average number of recurrent event per subject:  2.1 
#> Proportion of subjects with a terminal event:   0.2 
#> 
#> 
head(dat$data)
#>   id    t.start     t.stop event status        x1         x2        x3
#> 1  1 0.00000000 0.04635022     1      0 0.6273156 -0.4036852 0.6022021
#> 2  1 0.04635022 0.04805637     1      0 0.6273156 -0.4036852 0.6022021
#> 3  1 0.04805637 0.55984159     1      0 0.6273156 -0.4036852 0.6022021
#> 4  1 0.55984159 1.32060851     1      0 0.6273156 -0.4036852 0.6022021
#> 5  1 1.32060851 1.60517574     1      0 0.6273156 -0.4036852 0.6022021
#> 6  1 1.60517574 1.89312387     0      0 0.6273156 -0.4036852 0.6022021
#>          x4        x5        x6         x7         x8        x9       x10
#> 1 0.2667309 0.8732985 0.9716343 -0.5239045 -0.3774046 0.4123239 0.6563171
#> 2 0.2667309 0.8732985 0.9716343 -0.5239045 -0.3774046 0.4123239 0.6563171
#> 3 0.2667309 0.8732985 0.9716343 -0.5239045 -0.3774046 0.4123239 0.6563171
#> 4 0.2667309 0.8732985 0.9716343 -0.5239045 -0.3774046 0.4123239 0.6563171
#> 5 0.2667309 0.8732985 0.9716343 -0.5239045 -0.3774046 0.4123239 0.6563171
#> 6 0.2667309 0.8732985 0.9716343 -0.5239045 -0.3774046 0.4123239 0.6563171
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
# }
```
