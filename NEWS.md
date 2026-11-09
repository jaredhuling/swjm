# swjm 0.1.3

* The full joint frailty model (`estimate_frailty = TRUE`) now estimates the
  frailty variance. It was previously held at 1. The variance is updated at
  each path iteration by the closed-form solution of the estimating equation
  of Kalbfleisch et al. (2013, eq. 15), with the frailty weights held at their
  current values, and is returned as `frailty_var` by `stagewise_fit()` and
  `cv_stagewise()`.
* The frailty weights are predictable, `1 / (1 + theta * Lambda^D(t-))`, and
  the death baseline is solved exactly in time order given those weights.
* Fixed the weighted recurrence risk sets, in which a subject's weight reset
  to 1 after each recurrence until the next death time.
* `cv_stagewise()` scores held-out folds of the full model with the
  frailty-weighted estimating equations, using each fold's coefficients and
  frailty variance.
* `generate_data()` (JFM, `cov_type = "fixed"` or `"external"`) no longer
  drops a censored subject's last recurrence. New argument `cov_scale` sets
  the covariate standard deviation.
* The frailty-weighted passes track only the scalar risk-set sum and form the
  score with one matrix product, so their cost per death time is O(n) rather
  than O(n p).

* Fixed the adaptive step-size rule. The step was previously derived from
  `count_digits()`, which counted decimal places rather than magnitude. Two
  consequences: a dual norm of 1.27 was scored like 0.5 and received a step ten
  times smaller than intended, so paths starting above lambda = 1 advanced too
  slowly to traverse the path within `max_iter` (most visibly for the JFM, where
  coefficients were left roughly eight times too small); and once the dual norm
  fell below 1e-4 it printed in scientific notation, `count_digits()` returned
  0, and the step jumped back to 0.1, preventing the path from descending below
  lambda = 1e-4 at any iteration count.
* The step now starts at `eps` at the top of the path and is divided by ten each
  time the dual norm falls a further decade below its initial value. The rule is
  scale invariant and penalty independent.
* `eps` is now used. It was previously discarded by the adaptive branch for
  both model types and all three penalties, so the documented step-size argument
  had no effect.
* Removed `count_digits()`, which is no longer used.
* `stagewise_fit()` and `cv_stagewise()` gain `lambda_min_ratio`, which stops the
  path once the dual norm falls to that fraction of its value at the top of the
  path. With the step-size fix in place the JFM path reaches its solution in the
  first fifth of the iterations and then sits static: roughly 75% of the default
  5000 iterations were contributing under 1% of the total coefficient movement,
  while still inflating the lambda grid that cross-validation has to search.
  The default follows `glmnet` (0.01 when n > p, 1e-4 otherwise) and cuts JFM
  fit-plus-CV time by about 57% with coefficient estimates unchanged. The JSCM
  path is still moving at iteration 5000 and is left essentially untouched.
  Set `lambda_min_ratio = 0` for the previous behaviour.

# swjm 0.1.2

* Corrected the stagewise update directions under the scaled penalties. The
  argmax of the update objective over the scaled-penalty unit ball carries a
  factor 1/xi on every death-model coordinate (scaled lasso and mixed-sign
  cooperative branches: e_j sgn(g)/xi; group and same-sign cooperative
  branches: (g_a, g_b/xi^2)/||(g_a, g_b/xi)||_2). Earlier versions used
  unit-norm directions in rescaled gradient coordinates, which under-stepped
  the death-model block by the rescaling factor at every iteration and
  systematically over-shrank death-model coefficients.
* The gradient-rescaling factor is now frozen at its initial value, so each
  fit emulates a single fixed scaled penalty and the solution paths are well
  behaved.
