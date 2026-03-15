test_that("generate_data_jfm produces correct format", {
  set.seed(1)
  dat <- generate_data_jfm(n = 20, p = 10, scenario = 1)

  expect_named(dat, c("data", "alpha_true", "beta_true"))
  expect_s3_class(dat$data, "data.frame")

  expect_true("id" %in% names(dat$data))
  expect_true("t.start" %in% names(dat$data))
  expect_true("t.stop" %in% names(dat$data))
  expect_true("event" %in% names(dat$data))
  expect_true("status" %in% names(dat$data))
  expect_true("x1" %in% names(dat$data))
  expect_true("x10" %in% names(dat$data))

  expect_equal(length(dat$alpha_true), 10)
  expect_equal(length(dat$beta_true), 10)
  expect_equal(length(unique(dat$data$id)), 20)
})

test_that("generate_data_jfm scenarios differ", {
  set.seed(1)
  d1 <- generate_data_jfm(n = 10, p = 10, scenario = 1)
  d2 <- generate_data_jfm(n = 10, p = 10, scenario = 2)
  expect_false(identical(d1$alpha_true, d2$alpha_true))
})

test_that("generate_data dispatches correctly", {
  set.seed(1)
  d_jfm <- generate_data(n = 10, p = 10, scenario = 1, model = "jfm")
  expect_named(d_jfm, c("data", "alpha_true", "beta_true"))
})

test_that("generate_data_jscm produces correct format", {
  skip_if_not_installed("reReg")
  set.seed(1)
  dat <- generate_data_jscm(n = 20, p = 10, scenario = 1)

  expect_named(dat, c("data", "alpha_true", "beta_true"))
  expect_true("id" %in% names(dat$data))
  expect_true("x1" %in% names(dat$data))
  expect_equal(length(dat$alpha_true), 10)
  expect_equal(length(dat$beta_true), 10)
})
