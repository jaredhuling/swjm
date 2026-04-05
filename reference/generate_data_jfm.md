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
#>   id   t.start    t.stop event status          x1         x2         x3
#> 1  1 0.0000000 0.7548797     1      0 -1.25150957  0.5284880 -1.2476163
#> 2  1 0.7548797 0.7864475     1      0  0.02544607 -0.6784517 -0.3720642
#> 3  1 0.7864475 1.4852963     1      0 -0.38487961  1.0167902 -0.5236225
#> 4  1 1.4852963 1.6023504     0      1 -0.10524807  0.4496792  0.4996052
#> 5  2 0.0000000 0.6442640     1      0 -0.30080952  0.3112786 -0.4458649
#> 6  2 0.6442640 0.6723606     1      0  1.02350289 -0.1697170 -0.8738294
#>            x4          x5          x6         x7         x8          x9
#> 1 -0.04165134 -1.05473729 -0.09817368  0.2614318  0.9101125  1.70788034
#> 2  1.21136061 -0.25107658 -0.54145726 -0.0933376 -0.5413840  0.07153290
#> 3  0.02745382  0.57984851  1.33528304  0.3586251  1.5526201  1.56769961
#> 4 -0.18600815  1.24476382 -2.63185225 -2.7710264 -1.0185838  0.67406362
#> 5 -0.53930173 -0.56814019  0.98709927 -2.0407828 -0.1063646  1.72622018
#> 6  0.23936987 -0.02458833  0.85857453  1.3034761 -0.3886719 -0.02542339
#>          x10
#> 1 -0.1168501
#> 2 -1.8226220
#> 3  0.7985797
#> 4  1.1291482
#> 5  0.4177740
#> 6  0.2447490
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
```
