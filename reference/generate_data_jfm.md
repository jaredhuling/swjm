# Generate Simulated Data for the Joint Frailty Model (JFM)

Generates recurrent-event and terminal-event data under a Cox-type joint
frailty model. Ported from `Data_Generation_time_dependent_new()`.

## Usage

``` r
generate_data_jfm(
  n,
  p,
  scenario = 1,
  b = 6.5,
  lambda0_d = 0.041,
  lambda0_r = 1
)
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
  6.50).

- lambda0_d:

  Numeric. Baseline hazard rate for the terminal event (default 0.041).

- lambda0_r:

  Numeric. Baseline hazard rate for recurrent events (default 1).

## Value

A list with components:

- data:

  Data frame with columns `id`, `t.start`, `t.stop`, `event`, `status`,
  `x1`, ..., `xp`.

- alpha_true:

  True alpha (terminal) coefficients.

- beta_true:

  True beta (recurrence) coefficients.

## Details

Internally the simulation uses the Rondeau et al. (2007) convention
where `alpha` governs death and `beta` governs recurrence. The returned
`alpha_true` and `beta_true` are relabelled to match the package-wide
convention:

- `alpha_true`: recurrence (readmission) coefficients.

- `beta_true`: terminal (death) coefficients.

Within each subject the covariates are regenerated at each gap time,
yielding time-dependent covariates. Censoring times are `Uniform(1, b)`.

## Examples

``` r
dat <- generate_data_jfm(n = 30, p = 10, scenario = 1)
head(dat$data)
#>   id   t.start    t.stop event status          x1         x2          x3
#> 1  1 0.0000000 0.3841453     1      0  0.26535479  0.6190669  0.47108084
#> 2  1 0.3841453 2.0762212     1      0 -0.66046112 -0.1490251 -0.87548530
#> 3  1 2.0762212 2.1639101     1      0 -0.63081925 -0.9535960 -0.07998247
#> 4  1 2.1639101 2.2749376     1      0 -0.54249536 -2.0945782  0.37805006
#> 5  1 2.2749376 6.2244180     0      0 -0.05013821  2.2803904  1.37981056
#> 6  2 0.0000000 0.5320308     1      0  0.04422852  0.6270975 -0.94724955
#>           x4          x5          x6          x7         x8         x9
#> 1  0.4631442 -0.16016136 -0.03568547 -0.57468482  0.8453500 -0.5817661
#> 2  1.6990016  0.06819026  0.53318758 -0.08094899  2.4467814  0.3838990
#> 3  1.0626959  1.35105684  0.23329022  0.33389879 -0.7647703  0.7693845
#> 4 -0.7059342 -0.58710453  0.10008636  1.20206636  0.6598440 -0.6267454
#> 5 -0.5736360 -1.41683743 -1.24490746  0.90008248  1.9894303 -2.7869753
#> 6 -1.1228269 -0.02032365  1.22502554  0.08305898 -1.8859152  1.0569271
#>          x10
#> 1 -0.5023571
#> 2 -1.0521716
#> 3 -0.2150130
#> 4 -0.5297761
#> 5 -0.0231353
#> 6  0.6129122
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
```
