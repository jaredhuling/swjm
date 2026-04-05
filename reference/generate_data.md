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
#>   id    t.start     t.stop event status         x1          x2          x3
#> 1  1 0.00000000 5.08425244     0      0 -1.2209602 -1.05669143 -0.99689055
#> 2  2 0.00000000 0.03480276     1      0 -0.1742999 -0.31637148 -0.49300189
#> 3  2 0.03480276 3.23552656     0      0 -2.0679667 -0.08573976  0.03228655
#> 4  3 0.00000000 0.10750988     0      1 -0.4694870  0.11801278  1.24744155
#> 5  4 0.00000000 0.35064330     1      0 -0.4466444 -0.09540038  1.25076223
#> 6  4 0.35064330 0.81418576     1      0  0.3921569  0.44287977  0.65029280
#>           x4         x5          x6         x7         x8         x9
#> 1 -0.6726051  0.6652684  2.23607768 -0.6017982 -0.2916177 -0.8409749
#> 2  0.1907785 -1.1037450 -0.10641331  1.4589927 -0.3074743  0.7571814
#> 3  0.2896864  0.8089381  0.32626247 -1.7038575  0.5255046 -0.8030338
#> 4 -0.6883973 -1.3835696 -0.04293242 -0.2092342 -1.4669935  2.2357784
#> 5  0.5174422 -0.8772400 -0.61067715  1.6350373 -0.9457495  0.8412654
#> 6  0.2483775  1.0043266 -0.92091419 -0.3537209  0.6011269 -0.1243730
#>           x10
#> 1  2.03169952
#> 2 -1.71082539
#> 3 -0.23765572
#> 4 -0.62514496
#> 5  0.04248508
#> 6  1.08907050

# JSCM data with 30 subjects and 10 covariates
dat_jscm <- generate_data(n = 30, p = 10, scenario = 1, model = "jscm")
#> Call: 
#> reReg::simGSC(n = n, summary = TRUE, para = para, xmat = X, censoring = C, 
#>     frailty = gamma, tau = 60)
#> 
#> Summary:
#> Sample size:                                    30 
#> Number of recurrent event observed:             36 
#> Average number of recurrent event per subject:  1.2 
#> Proportion of subjects with a terminal event:   0.233 
#> 
#> 
head(dat_jscm$data)
#>   id    t.start     t.stop event status         x1         x2         x3
#> 1  1 0.00000000 0.25837665     1      0  0.8887092 -0.5173527 -0.2413731
#> 2  1 0.25837665 2.80992943     1      0  0.8887092 -0.5173527 -0.2413731
#> 3  1 2.80992943 3.25463114     0      0  0.8887092 -0.5173527 -0.2413731
#> 4  2 0.00000000 1.54721023     0      0 -0.8962707  0.3564506 -0.2278025
#> 5  3 0.00000000 0.02497537     1      0 -0.2942717 -0.6263212 -0.4358789
#> 6  3 0.02497537 0.11012372     1      0 -0.2942717 -0.6263212 -0.4358789
#>           x4         x5        x6         x7          x8         x9        x10
#> 1 -0.1931633  0.8857679 0.2149760 0.87070651  0.30583614 -0.7820145  0.1663135
#> 2 -0.1931633  0.8857679 0.2149760 0.87070651  0.30583614 -0.7820145  0.1663135
#> 3 -0.1931633  0.8857679 0.2149760 0.87070651  0.30583614 -0.7820145  0.1663135
#> 4  0.4469804 -0.3620662 0.7094065 0.33687002 -0.67381301 -0.8212125 -0.5228316
#> 5  0.9215679 -0.4508393 0.1236785 0.03369255  0.07984331  0.3145124 -0.2548009
#> 6  0.9215679 -0.4508393 0.1236785 0.03369255  0.07984331  0.3145124 -0.2548009
```
