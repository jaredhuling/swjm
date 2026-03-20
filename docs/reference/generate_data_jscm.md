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
#> Number of recurrent event observed:             39 
#> Average number of recurrent event per subject:  1.3 
#> Proportion of subjects with a terminal event:   0.167 
#> 
#> 
head(dat$data)
#>   id   t.start    t.stop event status         x1           x2         x3
#> 1  1 0.0000000 2.2434501     0      0  0.5908649 -0.128562628 -0.4892124
#> 2  2 0.0000000 3.1711364     0      0 -0.7642095 -0.164923138  0.2163843
#> 3  3 0.0000000 2.7547372     0      0  0.4209210  0.823433789  0.2011080
#> 4  4 0.0000000 0.9994934     1      0  0.7849151  0.503924026  0.4963014
#> 5  4 0.9994934 3.7459081     0      0  0.7849151  0.503924026  0.4963014
#> 6  5 0.0000000 0.1799844     1      0  0.4145064 -0.004255208 -0.8271047
#>           x4          x5         x6         x7         x8          x9
#> 1  0.4407820 -0.58694603 -0.4459489 -0.8577905  0.3268330 -0.36961792
#> 2  0.2424525 -0.07861706  0.4945090 -0.2397139 -0.9787622 -0.71564185
#> 3 -0.5282536 -0.88805921 -0.7635521 -0.1673474  0.2561140 -0.01826686
#> 4 -0.9305625 -0.75849909 -0.4077495  0.9353851 -0.6032338  0.30082858
#> 5 -0.9305625 -0.75849909 -0.4077495  0.9353851 -0.6032338  0.30082858
#> 6  0.2489779 -0.72096662  0.7699978  0.3252360 -0.6647410 -0.18642561
#>          x10
#> 1  0.9377569
#> 2  0.2787154
#> 3  0.2416653
#> 4  0.1977258
#> 5  0.1977258
#> 6 -0.9389974
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
# }
```
