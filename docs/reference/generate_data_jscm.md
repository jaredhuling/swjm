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
#> Proportion of subjects with a terminal event:   0.2 
#> 
#> 
head(dat$data)
#>   id   t.start    t.stop event status        x1         x2         x3        x4
#> 1  1 0.0000000 0.3222933     1      0 0.1313358 -0.4590607 -0.2478339 0.4998452
#> 2  1 0.3222933 0.5766185     1      0 0.1313358 -0.4590607 -0.2478339 0.4998452
#> 3  1 0.5766185 0.9428153     1      0 0.1313358 -0.4590607 -0.2478339 0.4998452
#> 4  1 0.9428153 1.2662159     1      0 0.1313358 -0.4590607 -0.2478339 0.4998452
#> 5  1 1.2662159 1.4086916     1      0 0.1313358 -0.4590607 -0.2478339 0.4998452
#> 6  1 1.4086916 2.2472515     1      0 0.1313358 -0.4590607 -0.2478339 0.4998452
#>          x5         x6         x7        x8        x9        x10
#> 1 0.8365068 -0.3354045 -0.1285257 0.5892879 0.6958139 -0.3235715
#> 2 0.8365068 -0.3354045 -0.1285257 0.5892879 0.6958139 -0.3235715
#> 3 0.8365068 -0.3354045 -0.1285257 0.5892879 0.6958139 -0.3235715
#> 4 0.8365068 -0.3354045 -0.1285257 0.5892879 0.6958139 -0.3235715
#> 5 0.8365068 -0.3354045 -0.1285257 0.5892879 0.6958139 -0.3235715
#> 6 0.8365068 -0.3354045 -0.1285257 0.5892879 0.6958139 -0.3235715
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
# }
```
