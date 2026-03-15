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
#> Number of recurrent event observed:             55 
#> Average number of recurrent event per subject:  1.833 
#> Proportion of subjects with a terminal event:   0.233 
#> 
#> 
head(dat$data)
#>   id   t.start    t.stop event status        x1        x2        x3        x4
#> 1  1 0.0000000 0.1073653     1      0 0.8375633 0.4189464 0.6367315 0.9597294
#> 2  1 0.1073653 0.2176367     1      0 0.8375633 0.4189464 0.6367315 0.9597294
#> 3  1 0.2176367 0.3734865     1      0 0.8375633 0.4189464 0.6367315 0.9597294
#> 4  1 0.3734865 0.9303313     1      0 0.8375633 0.4189464 0.6367315 0.9597294
#> 5  1 0.9303313 1.2339681     1      0 0.8375633 0.4189464 0.6367315 0.9597294
#> 6  1 1.2339681 2.3983258     0      0 0.8375633 0.4189464 0.6367315 0.9597294
#>           x5        x6         x7       x8         x9        x10
#> 1 -0.9885748 0.3998248 -0.8045619 0.270392 -0.3931547 -0.6490823
#> 2 -0.9885748 0.3998248 -0.8045619 0.270392 -0.3931547 -0.6490823
#> 3 -0.9885748 0.3998248 -0.8045619 0.270392 -0.3931547 -0.6490823
#> 4 -0.9885748 0.3998248 -0.8045619 0.270392 -0.3931547 -0.6490823
#> 5 -0.9885748 0.3998248 -0.8045619 0.270392 -0.3931547 -0.6490823
#> 6 -0.9885748 0.3998248 -0.8045619 0.270392 -0.3931547 -0.6490823
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
# }
```
