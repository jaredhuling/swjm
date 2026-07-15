# Changelog

## swjm 0.1.2

- Corrected the stagewise update directions under the scaled penalties.
  The argmax of the update objective over the scaled-penalty unit ball
  carries a factor 1/xi on every death-model coordinate (scaled lasso
  and mixed-sign cooperative branches: e_j sgn(g)/xi; group and
  same-sign cooperative branches: (g_a, g_b/xi^2)/\|\|(g_a,
  g_b/xi)\|\|\_2). Earlier versions used unit-norm directions in
  rescaled gradient coordinates, which under-stepped the death-model
  block by the rescaling factor at every iteration and systematically
  over-shrank death-model coefficients.
- The gradient-rescaling factor is now frozen at its initial value, so
  each fit emulates a single fixed scaled penalty and the solution paths
  are well behaved.
