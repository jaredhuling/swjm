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
#>   id    t.start     t.stop event status          x1         x2         x3
#> 1  1 0.00000000 0.08852674     1      0 -0.04283746 -1.2815593  0.9675822
#> 2  1 0.08852674 5.88543210     0      0 -1.45677852  0.2343405 -0.9515281
#> 3  2 0.00000000 0.02702974     1      0 -0.05814550  1.4905300  0.2675769
#> 4  2 0.02702974 1.96994382     0      0 -2.20710605  0.8953441 -0.5710303
#> 5  3 0.00000000 4.66633902     1      0  0.85504041  0.9297268 -0.5114203
#> 6  3 4.66633902 4.77357336     1      0 -0.11534663  0.4029612 -1.3295532
#>           x4         x5         x6         x7          x8         x9        x10
#> 1  1.0290489 -2.1664253 -0.3033482  0.1792680  1.42739582 -0.6199466 -0.1318424
#> 2  1.8984394 -1.0427705 -1.3337463  1.7647248  0.50606343 -0.6552922  0.6675088
#> 3 -0.3105904 -0.1514627  0.0889827 -1.4503736  0.84298116  1.3602523 -1.4819821
#> 4  0.8030739 -0.1511729 -0.1165839 -0.3122022 -0.08643189 -0.3113840 -0.5994553
#> 5  0.5756067 -1.1797607  0.8739227  0.7526550  0.63195295  0.4890517  0.4934662
#> 6  0.8500323  1.1056819  0.4041419 -0.1231415  0.44884772  2.0681537 -1.2495368
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
```
