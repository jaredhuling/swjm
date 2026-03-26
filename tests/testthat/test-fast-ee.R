# Tests for fast O(n_pseudo) JFM EE evaluation
# Verifies that jfm_s0s1_fast_cpp produces identical results to jfm_s0s1_cpp

test_that("fast S0t/S1t matches old implementation for death sub-model", {
  set.seed(42)
  dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
  dc <- swjm:::extract_data_components(dat$data)
  p <- 10L

  A <- swjm:::jfm_wt_death(1, numeric(p), dc$t.start, dc$I, dc$Z, dc$td,
                             rep(0.2, length(dc$td)), dc$td.id)
  pseudo <- A$pseudo_entries
  Z_pseudo <- pseudo[, 3:ncol(pseudo), drop = FALSE]
  td_sorted <- sort(dc$td)
  wt_m <- matrix(1, dc$n, length(dc$td))

  exit_times <- swjm:::jfm_compute_exit_times(pseudo, dc$Y)
  is_last <- as.integer(swjm:::jfm_compute_is_last(pseudo))
  tl <- swjm:::jfm_precompute_timeline(pseudo[, 1], exit_times, is_last,
                                         td_sorted, 0L)

  for (trial in 1:3) {
    coef <- rnorm(p, sd = 0.3)
    old <- swjm:::jfm_s0s1_cpp(dc$Y, wt_m, td_sorted,
                                 A$index_death_matrix, Z_pseudo, coef)
    new <- swjm:::jfm_s0s1_fast_cpp(tl$type, tl$idx, tl$size,
                                      Z_pseudo, coef, length(dc$td), dc$n)
    expect_equal(old$S0t, new$S0t, tolerance = 1e-12)
    expect_equal(old$S1t, new$S1t, tolerance = 1e-12)
  }
})


test_that("fast S0t/S1t matches old implementation for recurrent sub-model", {
  set.seed(42)
  dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
  dc <- swjm:::extract_data_components(dat$data)
  p <- 10L

  A <- swjm:::jfm_wt_death(1, numeric(p), dc$t.start, dc$I, dc$Z, dc$td,
                             rep(0.2, length(dc$td)), dc$td.id)
  pseudo <- A$pseudo_entries
  Z_pseudo <- pseudo[, 3:ncol(pseudo), drop = FALSE]
  tr_sorted <- sort(dc$tr)
  wt_r <- matrix(1, dc$n, length(dc$tr))

  B <- swjm:::jfm_r2i_integral(dc$t.start, dc$I, dc$Z, numeric(p),
                                 dc$tr, diff(c(0, tr_sorted)), dc$tr.id)

  exit_times <- swjm:::jfm_compute_exit_times(pseudo, dc$Y)
  is_last <- as.integer(swjm:::jfm_compute_is_last(pseudo))
  tl <- swjm:::jfm_precompute_timeline(pseudo[, 1], exit_times, is_last,
                                         tr_sorted, 1L)

  for (trial in 1:3) {
    coef <- rnorm(p, sd = 0.3)
    old <- swjm:::jfm_s0s1_cpp(dc$Y, wt_r, tr_sorted,
                                 B$index_recurrent_matrix, Z_pseudo, coef)
    new <- swjm:::jfm_s0s1_fast_cpp(tl$type, tl$idx, tl$size,
                                      Z_pseudo, coef, length(dc$tr), dc$n)
    expect_equal(old$S0t, new$S0t, tolerance = 1e-12)
    expect_equal(old$S1t, new$S1t, tolerance = 1e-12)
  }
})


test_that("fast score matches old score", {
  set.seed(42)
  dat <- generate_data(n = 50, p = 10, scenario = 1, model = "jfm")
  dc <- swjm:::extract_data_components(dat$data)
  p <- 10L

  A <- swjm:::jfm_wt_death(1, numeric(p), dc$t.start, dc$I, dc$Z, dc$td,
                             rep(0.2, length(dc$td)), dc$td.id)
  pseudo <- A$pseudo_entries
  Z_pseudo <- pseudo[, 3:ncol(pseudo), drop = FALSE]
  td_sorted <- sort(dc$td)
  wt_m <- matrix(1, dc$n, length(dc$td))

  exit_times <- swjm:::jfm_compute_exit_times(pseudo, dc$Y)
  is_last <- as.integer(swjm:::jfm_compute_is_last(pseudo))
  tl <- swjm:::jfm_precompute_timeline(pseudo[, 1], exit_times, is_last,
                                         td_sorted, 0L)
  de_epi <- swjm:::jfm_event_pseudo_idx(A$index_death_matrix, A$td_id)

  coef <- rnorm(p, sd = 0.3)
  res <- swjm:::jfm_s0s1_fast_cpp(tl$type, tl$idx, tl$size,
                                    Z_pseudo, coef, length(dc$td), dc$n)

  old_score <- swjm:::jfm_score_cpp(A$index_death_matrix, A$td_id - 1L,
                                      Z_pseudo, res$S1t, res$S0t)
  new_score <- swjm:::jfm_score_fast_cpp(de_epi, Z_pseudo, res$S1t, res$S0t)
  expect_equal(drop(old_score), drop(new_score), tolerance = 1e-12)
})


test_that("stagewise_fit produces consistent results with fast EE", {
  set.seed(99)
  dat <- generate_data(n = 30, p = 10, scenario = 1, model = "jfm")
  fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
                       max_iter = 50)

  # Should complete without error and produce reasonable output

  expect_s3_class(fit, "swjm_path")
  expect_equal(fit$k, 50L)
  expect_true(all(is.finite(fit$lambda)))
  expect_true(all(is.finite(fit$theta)))
})
