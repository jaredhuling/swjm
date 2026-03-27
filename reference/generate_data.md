# Generate Simulated Data for Joint Models

Unified interface that dispatches to model-specific data generation
functions for the joint frailty model (JFM) or joint scale-change model
(JSCM).

## Usage

``` r
generate_data(n, p, scenario = 1, model = c("jfm", "jscm"), ...)
```

## Arguments

- n:

  Integer. Number of subjects.

- p:

  Integer. Number of covariates (should be a multiple of 10 for
  scenarios 1–3).

- scenario:

  Integer. Scenario for true coefficient configuration (1, 2, 3, or
  other for a simple default).

- model:

  Character. Either `"jfm"` for the joint frailty model or `"jscm"` for
  the joint scale-change model.

- ...:

  Additional arguments passed to the model-specific function. For JFM:
  `b`, `lambda0_d`, `lambda0_r`, `gamma_frailty`, `cov_type`. For JSCM:
  `b`.

## Value

A list with components:

- data:

  A data frame in recurrent-event format with columns `id`, `t.start`,
  `t.stop`, `event`, `status`, and covariate columns `x1`, ..., `xp`.

- alpha_true:

  Numeric vector of true alpha coefficients.

- beta_true:

  Numeric vector of true beta coefficients.

## Examples

``` r
# JFM data with 30 subjects and 10 covariates
dat_jfm <- generate_data(n = 30, p = 10, scenario = 1, model = "jfm")
head(dat_jfm$data)
#>   id   t.start    t.stop event status       x1         x2         x3         x4
#> 1  1 0.0000000 0.0702400     1      0 1.678755 -0.2514572 -0.8552489  0.0293183
#> 2  1 0.0702400 0.9514857     1      0 0.814653 -0.4322850  1.1543770  1.5152614
#> 3  1 0.9514857 0.9942508     1      0 1.036150  0.2914693 -1.3300558  1.1996116
#> 4  1 0.9942508 1.0346228     1      0 0.860862 -0.7245456 -1.8403815  0.2428810
#> 5  1 1.0346228 1.0797941     1      0 1.831734 -0.9688322 -1.2200304 -0.2433400
#> 6  1 1.0797941 1.1003310     0      1 1.451714  1.3859527  0.9360747  1.0786122
#>           x5         x6         x7         x8         x9        x10
#> 1 -0.3520895 -1.3680452  0.4846134 -0.7003770 -0.1648789  0.4580152
#> 2  0.2178935  0.3404476  1.2930919  1.6066930 -0.4735040  0.3442692
#> 3 -0.1303352 -0.4299990 -0.9468287 -0.9769485  0.2561068  0.5028627
#> 4 -1.4133578 -1.7309186 -0.3731065  0.3566295  1.6890087 -0.2419267
#> 5  1.3764592  0.7078778 -2.0243059 -0.8033473 -0.1740919  0.5078751
#> 6 -0.5454287 -0.2392256  0.9703970  0.5309693  0.2837885 -0.5422261

# JSCM data with 30 subjects and 10 covariates
dat_jscm <- generate_data(n = 30, p = 10, scenario = 1, model = "jscm")
#> Call: 
#> reReg::simGSC(n = n, summary = TRUE, para = para, xmat = X, censoring = C, 
#>     frailty = gamma, tau = 60)
#> 
#> Summary:
#> Sample size:                                    30 
#> Number of recurrent event observed:             35 
#> Average number of recurrent event per subject:  1.167 
#> Proportion of subjects with a terminal event:   0.067 
#> 
#> 
head(dat_jscm$data)
#>   id  t.start     t.stop event status         x1          x2         x3
#> 1  1 0.000000 0.56248784     0      1 -0.5517991 -0.06476328  0.9378932
#> 2  2 0.000000 0.04136184     0      0  0.3327376  0.23642109  0.6357420
#> 3  3 0.000000 1.66010126     0      0  0.7909077  0.62182775 -0.3879670
#> 4  4 0.000000 0.56446304     1      0 -0.7087017  0.04359532  0.5828784
#> 5  4 0.564463 1.18669926     1      0 -0.7087017  0.04359532  0.5828784
#> 6  4 1.186699 2.30562160     0      0 -0.7087017  0.04359532  0.5828784
#>           x4         x5          x6         x7          x8          x9
#> 1 -0.9683988  0.7635903  0.36499811 -0.3731795  0.49704390 -0.91763139
#> 2  0.4920274 -0.4092843 -0.33273574 -0.7241345 -0.09685851 -0.03836687
#> 3  0.4294541 -0.9399628  0.07579597 -0.2790933 -0.34261811 -0.73352419
#> 4  0.3791268 -0.3704498 -0.48554485 -0.7020813 -0.19231589  0.15718414
#> 5  0.3791268 -0.3704498 -0.48554485 -0.7020813 -0.19231589  0.15718414
#> 6  0.3791268 -0.3704498 -0.48554485 -0.7020813 -0.19231589  0.15718414
#>           x10
#> 1  0.03391273
#> 2  0.80207595
#> 3  0.92720362
#> 4 -0.16906633
#> 5 -0.16906633
#> 6 -0.16906633
```
