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
#> Number of recurrent event observed:             50 
#> Average number of recurrent event per subject:  1.667 
#> Proportion of subjects with a terminal event:   0.133 
#> 
#> 
head(dat$data)
#>   id    t.start     t.stop event status        x1         x2         x3
#> 1  1 0.00000000 0.40241586     1      0 0.8729171 -0.3086353 -0.4350011
#> 2  1 0.40241586 3.01822092     0      0 0.8729171 -0.3086353 -0.4350011
#> 3  2 0.00000000 0.02667026     0      1 0.7378056 -0.1108706  0.9566610
#> 4  3 0.00000000 0.02808245     1      0 0.9356901  0.4935714 -0.6166567
#> 5  3 0.02808245 0.11825516     1      0 0.9356901  0.4935714 -0.6166567
#> 6  3 0.11825516 0.20618963     0      0 0.9356901  0.4935714 -0.6166567
#>           x4         x5         x6          x7         x8          x9
#> 1  0.5105003 -0.9849186  0.9945803  0.26765015  0.9199894 -0.04861799
#> 2  0.5105003 -0.9849186  0.9945803  0.26765015  0.9199894 -0.04861799
#> 3 -0.6860576  0.5610226  0.5260299 -0.24362089  0.3972881 -0.41722783
#> 4  0.7038770  0.6784643 -0.5714827 -0.01726216 -0.7016461  0.19980463
#> 5  0.7038770  0.6784643 -0.5714827 -0.01726216 -0.7016461  0.19980463
#> 6  0.7038770  0.6784643 -0.5714827 -0.01726216 -0.7016461  0.19980463
#>           x10
#> 1 -0.09717754
#> 2 -0.09717754
#> 3 -0.89417586
#> 4 -0.41797094
#> 5 -0.41797094
#> 6 -0.41797094
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
# }
```
