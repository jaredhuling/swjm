test_that("cv_stagewise runs for JFM", {
  set.seed(42)
  dat <- generate_data(n = 40, p = 10, scenario = 1, model = "jfm")
  fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
                       max_iter = 50)

  cv_res <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
                         lambda_seq = fit$lambda, K = 3, max_iter = 50)

  expect_s3_class(cv_res, "swjm_cv")
  expect_true(cv_res$position_CF >= 1)
  expect_true(cv_res$position_CF <= length(cv_res$lambda_seq))
  expect_equal(length(cv_res$Scorenorm_crossfit), length(cv_res$lambda_seq))
  expect_equal(length(cv_res$Scorenorm_crossfit_re), length(cv_res$lambda_seq))
})

test_that("coef.swjm_cv extracts coefficients for JFM", {
  set.seed(42)
  dat <- generate_data(n = 40, p = 10, scenario = 1, model = "jfm")
  fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
                       max_iter = 50)

  cv_res <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
                         lambda_seq = fit$lambda, K = 3, max_iter = 50)

  theta <- coef(cv_res)
  expect_length(theta, 20)
  expect_true(is.numeric(theta))
})

test_that("print.swjm_cv works", {
  set.seed(42)
  dat <- generate_data(n = 40, p = 10, scenario = 1, model = "jfm")
  fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
                       max_iter = 50)

  cv_res <- cv_stagewise(dat$data, model = "jfm", penalty = "coop",
                         lambda_seq = fit$lambda, K = 3, max_iter = 50)

  expect_output(print(cv_res), "Cross-validation")
})

test_that("cv_stagewise runs for JSCM", {
  skip_if_not_installed("reReg")
  set.seed(42)
  dat <- generate_data(n = 40, p = 10, scenario = 1, model = "jscm")
  fit <- stagewise_fit(dat$data, model = "jscm", penalty = "coop",
                       max_iter = 30)

  cv_res <- cv_stagewise(dat$data, model = "jscm", penalty = "coop",
                         lambda_seq = fit$lambda, K = 3, max_iter = 30)

  expect_s3_class(cv_res, "swjm_cv")
  expect_true(cv_res$position_CF >= 1)
})
