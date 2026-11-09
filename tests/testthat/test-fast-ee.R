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


test_that("score-only frailty-weighted passes match the S1-based versions", {
  set.seed(20261108)
  dat <- generate_data(n = 80, p = 6, scenario = 1, model = "jfm",
                       gamma_frailty = 1, cov_type = "fixed", cov_scale = 0.2)
  dc <- swjm:::extract_data_components(dat$data)
  n <- dc$n; p <- 6L
  psd <- swjm:::jfm_build_pseudo_cpp(dc$t.start, dc$I, do.call(rbind, dc$Z),
                                     dc$Y, dc$td, dc$td.id, dc$tr, dc$tr.id,
                                     numeric(p), rep(0.2, length(dc$td)), 1)
  Zp <- psd$Z_pseudo
  tl_d <- swjm:::jfm_precompute_timeline(psd$entry_times, psd$exit_times,
                                         psd$is_last, psd$td_sorted, 0L)
  tl_r <- swjm:::jfm_precompute_timeline_wt(psd$entry_times, psd$exit_times,
                                            psd$is_last, psd$tr_sorted,
                                            psd$td_sorted, 1L)
  for (trial in 1:3) {
    el_d <- exp(drop(Zp %*% rnorm(p, sd = 0.3)))
    el_r <- exp(drop(Zp %*% rnorm(p, sd = 0.3)))
    th <- runif(n, 0.2, 2)
    fw_old <- swjm:::jfm_death_forward_cpp(tl_d$type, tl_d$idx, tl_d$size, Zp,
                                           psd$entry_subject, el_d, th,
                                           psd$td_sorted, n)
    fw_new <- swjm:::jfm_death_forward_score_cpp(tl_d$type, tl_d$idx, tl_d$size,
                                                 Zp, psd$entry_subject, el_d, th,
                                                 psd$td_sorted, n, psd$de_epi)
    expect_equal(fw_new$score,
                 drop(swjm:::jfm_score_fast_cpp(psd$de_epi, Zp, fw_old$S1t,
                                                fw_old$S0t)), tolerance = 1e-10)
    for (nm in c("S0t", "lambda0_d", "wt_pre", "wt_post"))
      expect_equal(fw_new[[nm]], fw_old[[nm]], tolerance = 1e-12)

    rr_old <- swjm:::jfm_s0s1_wt_lp_cpp(tl_r$type, tl_r$idx, tl_r$size, Zp,
                                        psd$entry_subject, el_r, fw_old$wt_post,
                                        length(dc$tr), n)
    rr_new <- swjm:::jfm_wt_score_lp_cpp(tl_r$type, tl_r$idx, tl_r$size, Zp,
                                         psd$entry_subject, el_r, fw_old$wt_post,
                                         length(dc$tr), n, psd$re_epi)
    expect_equal(rr_new$S0t, rr_old$S0t, tolerance = 1e-12)
    expect_equal(rr_new$score,
                 drop(swjm:::jfm_score_fast_cpp(psd$re_epi, Zp, rr_old$S1t,
                                                rr_old$S0t)), tolerance = 1e-10)
  }
})
