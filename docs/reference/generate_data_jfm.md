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
#>   id    t.start     t.stop event status         x1         x2         x3
#> 1  1 0.00000000 2.63870834     0      0 -0.7868053  1.7628243 -0.3161494
#> 2  2 0.00000000 4.48365849     0      0 -0.4530235  1.7880689  0.9152691
#> 3  3 0.00000000 0.03537555     1      0  0.6526741 -0.8326049 -0.5655195
#> 4  3 0.03537555 3.77251799     0      0  0.5374563  2.0874372  0.6391305
#> 5  4 0.00000000 0.93162224     1      0  0.8312447  0.3032645  1.4509583
#> 6  4 0.93162224 1.39465120     0      0  0.5236375  1.6375199  1.0096613
#>          x4          x5          x6         x7          x8          x9
#> 1 0.6698138  1.85595228 -0.43665065  1.6834879 -0.22829121  0.83479832
#> 2 1.7081569  0.17981169 -1.10235307 -1.5030560  0.01791467 -0.05772346
#> 3 0.6834522 -0.79358592  0.01365963  0.5286262 -1.59398976  0.60226905
#> 4 0.7616294 -0.99666085  1.70989586  1.4109858  1.14954578 -0.75740247
#> 5 0.7540735 -0.38826349  1.78433042  0.9936069  0.81263033 -0.17099835
#> 6 0.3659726  0.02377992 -0.91268873  0.1207479  0.21850265 -1.06359971
#>          x10
#> 1 -0.6273759
#> 2 -0.7095120
#> 3 -0.4233526
#> 4  0.5821299
#> 5  0.1998907
#> 6  1.2784384
dat$alpha_true
#>  [1]  1.1 -1.1  0.1 -0.1  0.0  0.0  0.0  0.0  1.0 -1.0
dat$beta_true
#>  [1]  0.1 -0.1  1.1 -1.1  0.0  0.0  0.0  0.0  1.0 -1.0
```
