# Plot Predicted Survival Curves and Predictor Contributions

Produces a figure for a `"swjm_pred"` object. The layout depends on the
model:

## Usage

``` r
# S3 method for class 'swjm_pred'
plot(
  x,
  which_subject = 1L,
  which_process = c("both", "readmission", "death"),
  threshold = 0,
  ...
)
```

## Arguments

- x:

  An object of class `"swjm_pred"`.

- which_subject:

  Integer. Index of the subject to highlight (default `1`).

- which_process:

  Character. Which sub-model(s) to plot: `"both"` (default),
  `"readmission"`, or `"death"`.

- threshold:

  Non-negative numeric. Only predictors whose absolute subject-specific
  contribution exceeds this value are shown in the bar chart. The
  default `0` suppresses variables with exactly zero contribution (i.e.,
  unselected variables whose coefficient is zero). Increase `threshold`
  to focus on the most impactful predictors.

- ...:

  Currently unused.

## Value

Invisibly returns `x`.

## Details

**JFM** (four panels):

1.  Readmission-free survival curves for all subjects, with the selected
    subject highlighted.

2.  Bar chart of readmission predictor contributions \\\hat\alpha_j
    z\_{ij}\\ (log-hazard scale).

3.  Mortality survival curves with the selected subject highlighted.

4.  Bar chart of death predictor contributions \\\hat\beta_j z\_{ij}\\.

**JSCM** (four panels):

1.  Recurrent-event survival curves (AFT model) with the selected
    subject highlighted. The panel title shows the total acceleration
    factor \\e^{\hat\alpha^\top z_i}\\.

2.  Bar chart of recurrence log time-acceleration contributions
    \\\hat\alpha_j z\_{ij}\\: positive = events sooner, negative =
    later.

3.  Mortality survival curves with the selected subject highlighted,
    total acceleration factor \\e^{\hat\beta^\top z_i}\\ in the title.

4.  Bar chart of terminal-event log time-acceleration contributions
    \\\hat\beta_j z\_{ij}\\.

In all cases bars represent *subject-specific* contributions
(coefficient \\\times\\ covariate value), not bare coefficients, so the
display correctly reflects how much each predictor shifts the log-hazard
(JFM) or log time-acceleration (JSCM) for this particular subject.

## Examples

``` r
# \donttest{
dat  <- generate_data(n = 50, p = 5, scenario = 1, model = "jfm")
#> Error in beta %*% z: non-conformable arguments
cv   <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
                     max_iter = 100)
#> Error: object 'dat' not found
newz <- matrix(rnorm(15), nrow = 3, ncol = 5)
pred <- predict(cv, newdata = newz)
#> Error: object 'cv' not found
plot(pred, which_subject = 2)
#> Error: object 'pred' not found
plot(pred, which_subject = 2, threshold = 0.05)
#> Error: object 'pred' not found
# }
```
