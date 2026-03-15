test_that("lambda_interp works with single lambda", {
  result <- lambda_interp(5.0, c(3.0, 4.0, 6.0))
  expect_equal(result$left, c(1L, 1L, 1L))
  expect_equal(result$right, c(1L, 1L, 1L))
  expect_equal(result$frac, c(1, 1, 1))
})

test_that("lambda_interp works with multiple lambdas", {
  lam <- c(10, 8, 5, 2, 1)
  result <- lambda_interp(lam, c(9, 6, 1))
  expect_length(result$left, 3)
  expect_length(result$right, 3)
  expect_length(result$frac, 3)
  expect_true(all(result$frac >= 0 & result$frac <= 1))
})

test_that("coop_norm gives L2 for same-sign pairs", {
  theta <- c(3, 4)
  expect_equal(coop_norm(theta, 1), 5)
})

test_that("coop_norm gives L1 for opposite-sign pairs", {
  theta <- c(3, -4)
  expect_equal(coop_norm(theta, 1), 7)
})

test_that("coop_norm handles multiple pairs", {
  # theta = c(3, -1, 4, 2), p = 2

  # pair 1: (3, 4) same sign -> sqrt(9+16) = 5

  # pair 2: (-1, 2) opposite sign -> 1 + 2 = 3
  theta <- c(3, -1, 4, 2)
  p <- 2
  result <- coop_norm(theta, p)
  expect_equal(result, 8)
})

test_that("count_digits works", {
  expect_equal(count_digits(0.0012), 3L)
  expect_equal(count_digits(0.5), 1L)
  expect_equal(count_digits(5), 0L)
})

test_that("create_folds creates K folds", {
  set.seed(1)
  folds <- create_folds(1:20, 5)
  expect_length(folds, 5)
  expect_equal(sort(unname(unlist(folds))), 1:20)
})

test_that("extract_data_components returns correct structure", {
  set.seed(1)
  dat <- generate_data_jfm(n = 20, p = 10, scenario = 1)
  dc <- extract_data_components(dat$data)
  expect_equal(dc$n, 20)
  expect_equal(dc$p, 10)
  expect_length(dc$Z, 20)
  expect_true(length(dc$td) > 0)
  expect_true(length(dc$tr) > 0)
  expect_equal(length(dc$Y), 20)
  expect_equal(length(dc$STATUS), 20)
})

test_that("extract_decreasing_indices works", {
  x <- c(10, 8, 9, 7, 6, 8, 5)
  idx <- extract_decreasing_indices(x)
  expect_true(all(diff(x[idx]) < 0) || length(idx) == 1)
})
