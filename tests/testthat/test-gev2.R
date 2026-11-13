#context("GEV functions 2")

## Example shape parameters

# Positive
xi1 <- 0.1
low <- -1 / xi1

# Zero
xi2 <- 0

# Negative
xi3 <- -1e-7
up <- -1 / xi3

## Example input vectors
# For testing pdf, log_pdf and cdf
xvec <- c(-Inf, 0, Inf, NA)
x1 <- c(low, xvec)
x2 <- xvec
x3 <- c(up, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("rgev, length of output", {
  expect_length(rgev(n = 1, shape = xi1), 1)
  expect_length(rgev(n = 100, shape = xi1), 100)
  expect_length(rgev(n = 0, shape = xi1), 0)
  expect_error(rgev(n = -2, shape = xi1))

  expect_length(rgev(n = 1, shape = xi2), 1)
  expect_length(rgev(n = 100, shape = xi2), 100)
  expect_length(rgev(n = 0, shape = xi2), 0)
  expect_error(rgev(n = -2, shape = xi2))

  expect_length(rgev(n = 1, shape = xi3), 1)
  expect_length(rgev(n = 100, shape = xi3), 100)
  expect_length(rgev(n = 0, shape = xi3), 0)
  expect_error(rgev(n = -2, shape = xi3))
})

test_that("dgev, values and length of output", {
  p <- pvec[2:4]
  qq <- qgev(p, shape = xi1)
  expect_equal(dgev(x1, shape = xi1), c(0, 0, exp(-1), 0, NA))
  expect_equal(dgev(qq, shape = xi1), (-log(p)) ^ (1 + xi1) * p)
  expect_length(dgev(seq_len(0), shape = xi1), 0)
  expect_length(dgev(seq_len(1), shape = xi1), 1)
  expect_length(dgev(seq_len(10), shape = xi1), 10)

  qq <- qgev(p, shape = xi2)
  expect_equal(dgev(x1, shape = xi2), c(0, 0, exp(-1), 0, NA))
  expect_equal(dgev(qq, shape = xi2), (-log(p)) ^ (1 + xi2) * p)
  expect_length(dgev(seq_len(0), shape = xi2), 0)
  expect_length(dgev(seq_len(1), shape = xi2), 1)
  expect_length(dgev(seq_len(10), shape = xi2), 10)

  qq <- qgev(p, shape = xi3)
  expect_equal(dgev(x1, shape = xi3), c(0, 0, exp(-1), 0, NA))
  expect_equal(dgev(qq, shape = xi3), (-log(p)) ^ (1 + xi3) * p)
  expect_length(dgev(seq_len(0), shape = xi3), 0)
  expect_length(dgev(seq_len(1), shape = xi3), 1)
  expect_length(dgev(seq_len(10), shape = xi3), 10)
})

test_that("dgev, allowing for NAs in x and all parameters (loc, scale, shape)", {
  ## Baseline, no missing values in x and any of the parameters
  expect_silent(baseline <- dgev(x = 0:4, loc = 0:4, scale = 1:5, shape = 0:4))

  ## Missing values in 'x'
  expect_silent(tmp <- dgev(x = c(NA, 1, NA, NA, 4), loc = 0:4, scale = 1:5, shape = 0:4))
  expect_identical(tmp, baseline * c(NA, 1, NA, NA, 1),
    info = "expected different return when 'x' includes missing values")

  ## Missing values in 'loc'
  expect_silent(tmp <- dgev(x = 0:4, loc = c(0, NA, 2, NA, NA), scale = 1:5, shape = 0:4))
  expect_identical(tmp, baseline * c(1, NA, 1, NA, NA),
    info = "expected different return when 'loc' includes missing values")

  ## Missing values in 'scale'
  expect_silent(tmp <- dgev(x = 0:4, loc = 0:4, scale = c(1, 2, NA, NA, 5), shape = 0:4))
  expect_identical(tmp, baseline * c(1, 1, NA, NA, 1),
    info = "expected different return when 'scale' includes missing values")

  ## Missing values in 'shape'
  expect_silent(tmp <- dgev(x = 0:4, loc = 0:4, scale = 1:5, shape = c(NA, 1, 2, 3, NA)))
  expect_identical(tmp, baseline * c(NA, 1, 1, 1, NA),
    info = "expected different return when 'shape' includes missing values")
})

test_that("dgev, allowing for NAs in x and all parameters (loc, scale, shape) w/ log = TRUE", {
  ## Baseline, no missing values in x and any of the parameters
  expect_silent(baseline <- dgev(x = 0:4, loc = 0:4, scale = 1:5, shape = 0:4, log = TRUE))

  ## Missing values in 'x'
  expect_silent(tmp <- dgev(x = c(NA, 1, NA, NA, 4), loc = 0:4, scale = 1:5, shape = 0:4, log = TRUE))
  expect_identical(tmp, baseline * c(NA, 1, NA, NA, 1),
    info = "expected different return when 'x' includes missing values")

  ## Missing values in 'loc'
  expect_silent(tmp <- dgev(x = 0:4, loc = c(0, NA, 2, NA, NA), scale = 1:5, shape = 0:4, log = TRUE))
  expect_identical(tmp, baseline * c(1, NA, 1, NA, NA),
    info = "expected different return when 'loc' includes missing values")

  ## Missing values in 'scale'
  expect_silent(tmp <- dgev(x = 0:4, loc = 0:4, scale = c(1, 2, NA, NA, 5), shape = 0:4, log = TRUE))
  expect_identical(tmp, baseline * c(1, 1, NA, NA, 1),
    info = "expected different return when 'scale' includes missing values")

  ## Missing values in 'shape'
  expect_silent(tmp <- dgev(x = 0:4, loc = 0:4, scale = 1:5, shape = c(NA, 1, 2, 3, NA), log = TRUE))
  expect_identical(tmp, baseline * c(NA, 1, 1, 1, NA),
    info = "expected different return when 'shape' includes missing values")
})

test_that("pgev, values and length of output", {
  expect_equal(pgev(x1, shape = xi1), c(0, 0, exp(-1), 1, NA))
  expect_length(pgev(seq_len(0), shape = xi1), 0)
  expect_length(pgev(seq_len(1), shape = xi1), 1)
  expect_length(pgev(seq_len(10), shape = xi1), 10)

  expect_equal(pgev(x2, shape = xi2), c(0, exp(-1), 1, NA))
  expect_length(pgev(seq_len(0), shape = xi2), 0)
  expect_length(pgev(seq_len(1), shape = xi2), 1)
  expect_length(pgev(seq_len(10), shape = xi2), 10)

  expect_equal(pgev(x3, shape = xi3), c(1, 0, exp(-1), 1, NA))
  expect_length(pgev(seq_len(0), shape = xi3), 0)
  expect_length(pgev(seq_len(1), shape = xi3), 1)
  expect_length(pgev(seq_len(10), shape = xi3), 10)
})

test_that("qgev, values and length of output", {
  q1 <- ((-log(pvec[2:4])) ^ (-xi1) - 1) / xi1
  expect_equal(qgev(pvec, shape = xi1), c(low, q1, Inf, NA))
  expect_length(qgev(seq_len(0), shape = xi1), 0)
  expect_length(qgev(c(0, 1), shape = xi1), 2)
  expect_length(qgev(seq_len(10) / 10, shape = xi1), 10)

  q2 <- -log(-log(pvec[2:4]))
  expect_equal(qgev(pvec, shape = xi2), c(-Inf, q2, Inf, NA))
  expect_length(qgev(seq_len(0), shape = xi2), 0)
  expect_length(qgev(c(0, 1), shape = xi2), 2)
  expect_length(qgev(seq_len(10) / 10, shape = xi2), 10)

  q3 <- ((-log(pvec[2:4])) ^ (-xi3) - 1) / xi3
  expect_equal(qgev(pvec, shape = xi3), c(-Inf, q3, up, NA))
  expect_length(qgev(seq_len(0), shape = xi3), 0)
  expect_length(qgev(c(0, 1), shape = xi3), 2)
  expect_length(qgev(seq_len(10) / 10, shape = xi3), 10)
})

test_that("pgev and qgev are consistent", {
  expect_equal(pgev(qgev(pvec, shape = xi1), shape = xi1), pvec)
  expect_equal(pgev(qgev(pvec, shape = xi2), shape = xi2), pvec)
  expect_equal(pgev(qgev(pvec, shape = xi3), shape = xi3), pvec)
})

test_that("qgev throws warning and returns NaN if p < 0 or p > 1", {
  # p < 0
  expect_warning(q <- qgev(-0.0001, shape = xi1), regexp = "NaNs produced")
  expect_identical(q, NaN)
  # p > 1
  expect_warning(q <- qgev(1.0001, shape = xi1), regexp = "NaNs produced")
  expect_identical(q, NaN)
  # exp(p) > 1
  expect_warning(q <- qgev(0.0001, shape = xi1, log.p = TRUE), regexp = "NaNs produced")
  expect_identical(q, NaN)

  p <- seq(-0.2, 1.2, by = 0.1)
  lidx <- p < 0 | p > 1
  expect_warning(q <- qgev(p, shape = xi1), regexp = "NaNs produced")
  expect_true(all(is.nan(q[lidx])))
  expect_true(all(!is.nan(q[!lidx])))

  p <- seq(-2, 1, by = 0.1)
  lidx <- p > 0
  expect_warning(q <- qgev(p, shape = xi1, log.p = TRUE), regexp = "NaNs produced")
  expect_true(all(is.nan(q[lidx])))
  expect_true(all(!is.nan(q[!lidx])))
})




