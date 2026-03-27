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
#>   id  t.start   t.stop event status          x1         x2          x3
#> 1  1 0.000000 1.766124     0      0  0.62181043  0.4134333 -0.09919527
#> 2  2 0.000000 1.237081     1      0  0.50181391 -1.2923657  0.81758778
#> 3  2 1.237081 1.253028     1      0  0.33158034  0.5087287 -0.49710559
#> 4  2 1.253028 1.627053     1      0  0.51382026  0.5027893 -0.73529158
#> 5  2 1.627053 1.710380     1      0 -0.01094071 -1.3926090  0.62494754
#> 6  2 1.710380 2.307124     1      0  1.00547316 -0.2973614 -0.23543516
#>           x4         x5          x6          x7         x8         x9
#> 1 -0.1960234  0.6211624 -0.04749709 -1.09567871  0.5978357 -0.8101280
#> 2 -0.7008723 -0.3608976 -0.68232175 -1.13663536 -0.5642234 -1.3490556
#> 3 -0.9443898 -1.1944597  1.35163231 -0.61388379 -0.6667017 -1.2925475
#> 4 -1.1085182  0.7170572  0.40585690 -1.49918152 -1.2235785  1.4253480
#> 5  0.6807196  2.2511400 -1.84212098 -0.08893203 -0.6541858 -0.3710897
#> 6 -0.2631708 -0.3204640  1.93065924 -2.18777451  0.7009740 -0.3763015
#>          x10
#> 1  0.8191176
#> 2 -0.6508801
#> 3  0.8227224
#> 4  0.5354840
#> 5  0.7321427
#> 6  0.3351810

# JSCM data with 30 subjects and 10 covariates
dat_jscm <- generate_data(n = 30, p = 10, scenario = 1, model = "jscm")
#> Call: 
#> reReg::simGSC(n = n, summary = TRUE, para = para, xmat = X, censoring = C, 
#>     frailty = gamma, tau = 60)
#> 
#> Summary:
#> Sample size:                                    30 
#> Number of recurrent event observed:             57 
#> Average number of recurrent event per subject:  1.9 
#> Proportion of subjects with a terminal event:   0.1 
#> 
#> 
head(dat_jscm$data)
#>   id   t.start    t.stop event status         x1         x2         x3       x4
#> 1  1 0.0000000 0.2164831     1      0 -0.1539435 -0.3572376 -0.3736188 0.304346
#> 2  1 0.2164831 0.6010202     1      0 -0.1539435 -0.3572376 -0.3736188 0.304346
#> 3  1 0.6010202 0.6150349     1      0 -0.1539435 -0.3572376 -0.3736188 0.304346
#> 4  1 0.6150349 1.0284986     1      0 -0.1539435 -0.3572376 -0.3736188 0.304346
#> 5  1 1.0284986 1.8690062     1      0 -0.1539435 -0.3572376 -0.3736188 0.304346
#> 6  1 1.8690062 2.2329656     1      0 -0.1539435 -0.3572376 -0.3736188 0.304346
#>           x5       x6        x7         x8        x9       x10
#> 1 -0.3084554 0.426761 0.4129865 -0.6660001 0.7307347 0.4135242
#> 2 -0.3084554 0.426761 0.4129865 -0.6660001 0.7307347 0.4135242
#> 3 -0.3084554 0.426761 0.4129865 -0.6660001 0.7307347 0.4135242
#> 4 -0.3084554 0.426761 0.4129865 -0.6660001 0.7307347 0.4135242
#> 5 -0.3084554 0.426761 0.4129865 -0.6660001 0.7307347 0.4135242
#> 6 -0.3084554 0.426761 0.4129865 -0.6660001 0.7307347 0.4135242
```
