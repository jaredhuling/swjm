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
  `b`, `lambda0_d`, `lambda0_r`. For JSCM: `b`.

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
#>   id    t.start     t.stop event status          x1         x2         x3
#> 1  1 0.00000000 3.97597637     1      0 -0.51485712 -0.8237882  0.3344152
#> 2  1 3.97597637 4.54024466     1      0 -0.12940880  0.4077955 -0.5839690
#> 3  1 4.54024466 5.29132630     0      1  0.17643641  1.1283074  0.4350017
#> 4  2 0.00000000 0.06902837     1      0 -0.06505926 -1.1300015  0.8250114
#> 5  2 0.06902837 0.48749959     1      0 -1.08769667 -1.3436343 -0.3595517
#> 6  2 0.48749959 2.60244354     0      0 -0.45335837  0.5229102  0.3534007
#>            x4         x5          x6         x7          x8          x9
#> 1 -0.09349076  0.3040621 -0.47650753 -0.2413117  0.82415583 -1.55564397
#> 2 -0.19344989 -0.2695467  0.07365823  0.3572361  0.55042842  0.03840179
#> 3  0.54883349  0.6474183  0.87846345  0.3507979  0.04987972  0.83574945
#> 4  0.59267509 -0.6603480  0.05987437 -1.0554508  0.98828970 -1.35048647
#> 5 -1.25745660 -1.1576458 -0.19059760 -0.5037913 -0.46928925  1.00168370
#> 6  0.22512154  1.6418438 -0.85181344 -1.8000535 -0.34388198 -2.13449444
#>           x10
#> 1  0.09350128
#> 2 -1.60957529
#> 3 -0.28172520
#> 4 -1.12441979
#> 5  1.27346352
#> 6 -0.13448241

# JSCM data with 30 subjects and 10 covariates
dat_jscm <- generate_data(n = 30, p = 10, scenario = 1, model = "jscm")
#> Call: 
#> reReg::simGSC(n = n, summary = TRUE, para = para, xmat = X, censoring = C, 
#>     frailty = gamma, tau = 60)
#> 
#> Summary:
#> Sample size:                                    30 
#> Number of recurrent event observed:             50 
#> Average number of recurrent event per subject:  1.667 
#> Proportion of subjects with a terminal event:   0.067 
#> 
#> 
head(dat_jscm$data)
#>   id t.start    t.stop event status         x1         x2          x3
#> 1  1       0 0.1388277     0      0 -0.2755367 -0.6828487 -0.23003906
#> 2  2       0 2.8749526     0      0 -0.7810648  0.3733265 -0.85125178
#> 3  3       0 0.6107033     0      0 -0.7184871  0.7565505 -0.26205813
#> 4  4       0 0.5986426     0      0  0.1984379 -0.6086902  0.02565082
#> 5  5       0 0.2644325     0      0 -0.4114530  0.5541339 -0.01368887
#> 6  6       0 3.7632460     0      0 -0.8819981 -0.1654656 -0.96580090
#>           x4          x5          x6          x7          x8         x9
#> 1  0.6158881 -0.25083423  0.03784226 -0.04373627  0.70898540 -0.7784821
#> 2 -0.1160849 -0.16373653  0.97508132 -0.69127673 -0.05017387  0.8925907
#> 3 -0.6389267  0.39283840 -0.35335505 -0.65802681  0.50926886 -0.2406596
#> 4  0.2589555  0.05974337  0.85775600 -0.25143003 -0.89820631  0.9168687
#> 5  0.1761063  0.89341738  0.24491763 -0.99651686  0.57337172 -0.8853154
#> 6 -0.2694233 -0.97560820  0.41909394 -0.35250114  0.75728586  0.7979397
#>           x10
#> 1  0.54767160
#> 2  0.96232522
#> 3  0.12276717
#> 4  0.94503274
#> 5 -0.25165767
#> 6  0.06020683
```
