# Changelog

## swjm 0.2.0

- Corrected the stagewise update directions under the scaled penalties.
  The argmax of the update objective over the scaled-penalty unit ball
  carries a factor 1/xi on every death-model coordinate (scaled lasso
  and mixed-sign cooperative branches: e_j sgn(g)/xi; group and
  same-sign cooperative branches: (g_a, g_b/xi^2)/\|\|(g_a,
  g_b/xi)\|\|\_2). The previous implementation used unit-norm directions
  in rescaled gradient coordinates, which under-stepped the death-model
  block by the rescaling factor at every iteration and systematically
  over-shrank death-model coefficients.
- New `direction` argument for
  [`stagewise_fit()`](http://jaredhuling.org/swjm/reference/stagewise_fit.md)
  and
  [`cv_stagewise()`](http://jaredhuling.org/swjm/reference/cv_stagewise.md)
  with variants `"corrected-fixed"` (default, recommended: exact
  directions with the rescaling factor frozen at initialization, giving
  well-behaved paths), `"corrected"` (exact directions with
  per-iteration recalibration), `"legacy"` (pre-correction behavior, for
  reproducibility), and `"corrected-capped"` (comparison only).
