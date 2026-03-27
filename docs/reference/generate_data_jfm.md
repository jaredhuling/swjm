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
  lambda0_r = 1,
  gamma_frailty = 0,
  cov_type = c("internal", "external", "fixed")
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

- gamma_frailty:

  Numeric. Frailty variance parameter. When positive, a subject-specific
  frailty \\Z_i \sim \text{Gamma}(1/\gamma, 1/\gamma)\\ is drawn for
  each subject and multiplies both hazard rates. When 0 (default), no
  frailty is used (\\Z_i = 1\\).

- cov_type:

  Character. How time-varying covariates are generated: `"internal"`
  (default) redraws covariates at each recurrent event; `"external"`
  changes covariates at predetermined Poisson times independent of the
  event process (Kalbfleisch-compatible); `"fixed"` draws one covariate
  vector per subject that never changes.

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
#>   id    t.start     t.stop event status         x1          x2         x3
#> 1  1 0.00000000 3.23340736     0      0 -1.0462034  0.17226001 -0.3070208
#> 2  2 0.00000000 0.08145319     1      0  0.9420595  0.46959178 -0.6621122
#> 3  2 0.08145319 0.24375363     1      0  0.4710808  0.46314421 -0.1601614
#> 4  2 0.24375363 0.39026596     1      0 -0.1490251 -0.87548530  1.6990016
#> 5  2 0.39026596 3.02040351     0      0 -0.9535960 -0.07998247  1.0626959
#> 6  3 0.00000000 5.42101568     0      0 -2.0945782  0.37805006 -0.7059342
#>            x4         x5          x6         x7         x8         x9
#> 1  0.90334561  0.9852487  1.64382848  0.8145597  0.8251132 -0.3457739
#> 2 -0.30947110  0.7945743 -0.91485096 -0.4655203  0.8196694  1.5210357
#> 3 -0.03568547 -0.5746848  0.84535002 -0.5817661 -0.5023571  1.0880763
#> 4  0.06819026  0.5331876 -0.08094899  2.4467814  0.3838990 -1.0521716
#> 5  1.35105684  0.2332902  0.33389879 -0.7647703  0.7693845 -0.2150130
#> 6 -0.58710453  0.1000864  1.20206636  0.6598440 -0.6267454 -0.5297761
#>          x10
#> 1  0.4111473
#> 2  0.2653548
#> 3 -0.3929530
#> 4 -1.7222338
#> 5  0.3833806
#> 6  0.8765032
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
```
