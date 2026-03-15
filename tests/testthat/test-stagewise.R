test_that("stagewise_fit runs for JFM coop", {
  set.seed(42)
  dat <- generate_data(n = 30, p = 10, scenario = 1, model = "jfm")
  fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
                       max_iter = 20)

  expect_s3_class(fit, "swjm_path")
  expect_equal(fit$model, "jfm")
  expect_equal(fit$penalty, "coop")
  expect_equal(fit$p, 10)
  expect_true(fit$k > 0)
  expect_equal(nrow(fit$theta), 20)
  expect_equal(ncol(fit$theta), fit$k + 1)
  expect_equal(length(fit$lambda), fit$k + 1)
})

test_that("stagewise_fit runs for JFM lasso", {
  set.seed(42)
  dat <- generate_data(n = 30, p = 10, scenario = 1, model = "jfm")
  fit <- stagewise_fit(dat$data, model = "jfm", penalty = "lasso",
                       max_iter = 20)
  expect_s3_class(fit, "swjm_path")
  expect_equal(fit$penalty, "lasso")
})

test_that("stagewise_fit runs for JFM group", {
  set.seed(42)
  dat <- generate_data(n = 30, p = 10, scenario = 1, model = "jfm")
  fit <- stagewise_fit(dat$data, model = "jfm", penalty = "group",
                       max_iter = 20)
  expect_s3_class(fit, "swjm_path")
  expect_equal(fit$penalty, "group")
})

test_that("stagewise_fit runs for JSCM coop", {
  skip_if_not_installed("reReg")
  set.seed(42)
  dat <- generate_data(n = 30, p = 10, scenario = 1, model = "jscm")
  fit <- stagewise_fit(dat$data, model = "jscm", penalty = "coop",
                       max_iter = 20)

  expect_s3_class(fit, "swjm_path")
  expect_equal(fit$model, "jscm")
  expect_equal(fit$p, 10)
  expect_true(fit$k > 0)
  expect_equal(nrow(fit$theta), 20)
})

test_that("stagewise_fit runs for JSCM lasso", {
  skip_if_not_installed("reReg")
  set.seed(42)
  dat <- generate_data(n = 30, p = 10, scenario = 1, model = "jscm")
  fit <- stagewise_fit(dat$data, model = "jscm", penalty = "lasso",
                       max_iter = 20)
  expect_s3_class(fit, "swjm_path")
})

test_that("print.swjm_path works", {
  set.seed(42)
  dat <- generate_data(n = 30, p = 10, scenario = 1, model = "jfm")
  fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
                       max_iter = 20)
  expect_output(print(fit), "Stagewise path")
})

test_that("theta starts at zero", {
  set.seed(42)
  dat <- generate_data(n = 30, p = 10, scenario = 1, model = "jfm")
  fit <- stagewise_fit(dat$data, model = "jfm", penalty = "coop",
                       max_iter = 20)
  expect_true(all(fit$theta[, 1] == 0))
})
