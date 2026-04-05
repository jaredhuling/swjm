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
#>   id    t.start     t.stop event status         x1         x2         x3
#> 1  1 0.00000000 0.47167890     1      0 -0.3112229 -0.1999889  1.6624539
#> 2  1 0.47167890 0.90309530     1      0 -0.1568004  1.1286446 -0.6893182
#> 3  1 0.90309530 3.59153677     0      0  1.6058849  0.9972221 -0.1124492
#> 4  2 0.00000000 0.07402463     1      0  0.7711456 -0.5861508  1.2759097
#> 5  2 0.07402463 0.50103401     1      0 -1.3574064  1.3463449 -0.5668461
#> 6  2 0.50103401 1.00240356     1      0 -1.1318343 -0.7923724 -0.2657469
#>           x4          x5         x6         x7         x8         x9        x10
#> 1 -0.3482454  0.86683165  1.0940562  0.9895212  0.9340026 -0.4137564 -0.2366623
#> 2 -0.0491129 -0.07360191  1.3229967 -0.7143528  2.2015747 -0.2526608  0.7971561
#> 3 -0.2780079 -0.18039598  1.2577244 -0.5554096 -0.9202618 -0.9579029  1.9574632
#> 4 -1.7784363 -1.36283784 -0.6789834 -1.2805472 -1.5696787 -1.9017783  0.8746232
#> 5 -0.3138270 -0.19934090 -1.4006258 -1.2317068 -0.3388905  0.3418990 -1.0224838
#> 6  0.6677651  0.95855022  1.9090526  0.7074147 -1.4047409  0.7940592 -0.4466444

# JSCM data with 30 subjects and 10 covariates
dat_jscm <- generate_data(n = 30, p = 10, scenario = 1, model = "jscm")
#> Call: 
#> reReg::simGSC(n = n, summary = TRUE, para = para, xmat = X, censoring = C, 
#>     frailty = gamma, tau = 60)
#> 
#> Summary:
#> Sample size:                                    30 
#> Number of recurrent event observed:             55 
#> Average number of recurrent event per subject:  1.833 
#> Proportion of subjects with a terminal event:   0.133 
#> 
#> 
head(dat_jscm$data)
#>   id     t.start      t.stop event status         x1         x2           x3
#> 1  1 0.000000000 0.002937668     1      0 0.83054436 -0.9898278 -0.009630135
#> 2  1 0.002937668 0.162268859     1      0 0.83054436 -0.9898278 -0.009630135
#> 3  1 0.162268859 0.504022744     1      0 0.83054436 -0.9898278 -0.009630135
#> 4  1 0.504022744 0.506660392     0      0 0.83054436 -0.9898278 -0.009630135
#> 5  2 0.000000000 0.134282175     1      0 0.04500048  0.3611698 -0.835212073
#> 6  2 0.134282175 0.469690210     1      0 0.04500048  0.3611698 -0.835212073
#>          x4         x5        x6        x7         x8         x9       x10
#> 1 0.3193714 -0.9651260 0.8754959 0.4624372 0.90183328 -0.2650918 0.3631836
#> 2 0.3193714 -0.9651260 0.8754959 0.4624372 0.90183328 -0.2650918 0.3631836
#> 3 0.3193714 -0.9651260 0.8754959 0.4624372 0.90183328 -0.2650918 0.3631836
#> 4 0.3193714 -0.9651260 0.8754959 0.4624372 0.90183328 -0.2650918 0.3631836
#> 5 0.5685206 -0.8585753 0.2556029 0.8921280 0.07164997  0.8693781 0.1886130
#> 6 0.5685206 -0.8585753 0.2556029 0.8921280 0.07164997  0.8693781 0.1886130
```
