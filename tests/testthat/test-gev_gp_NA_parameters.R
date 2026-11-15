# This stems from the first part of Issue 11
#   https://github.com/paulnorthrop/revdbayes/issues/11
# and the associated pull request 10
#   https://github.com/paulnorthrop/revdbayes/pull/10/changes
# It affects `dgev()`, `pgev()`, `qgev()`, `dgp()`, `pgp()` and `ggp()`

# GEV distribution functions

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

# GP distribution functions

test_that("dgp, allowing for NAs in x and all parameters (loc, scale, shape)", {
  ## Baseline, no missing values in x and any of the parameters
  expect_silent(baseline <- dgp(x = 0:4, loc = 0:4, scale = 1:5, shape = 0:4))

  ## Missing values in 'x'
  expect_silent(tmp <- dgp(x = c(NA, 1, NA, NA, 4), loc = 0:4, scale = 1:5, shape = 0:4))
  expect_identical(tmp, baseline * c(NA, 1, NA, NA, 1),
                   info = "expected different return when 'x' includes missing values")

  ## Missing values in 'loc'
  expect_silent(tmp <- dgp(x = 0:4, loc = c(0, NA, 2, NA, NA), scale = 1:5, shape = 0:4))
  expect_identical(tmp, baseline * c(1, NA, 1, NA, NA),
                   info = "expected different return when 'loc' includes missing values")

  ## Missing values in 'scale'
  expect_silent(tmp <- dgp(x = 0:4, loc = 0:4, scale = c(1, 2, NA, NA, 5), shape = 0:4))
  expect_identical(tmp, baseline * c(1, 1, NA, NA, 1),
                   info = "expected different return when 'scale' includes missing values")

  ## Missing values in 'shape'
  expect_silent(tmp <- dgp(x = 0:4, loc = 0:4, scale = 1:5, shape = c(NA, 1, 2, 3, NA)))
  expect_identical(tmp, baseline * c(NA, 1, 1, 1, NA),
                   info = "expected different return when 'shape' includes missing values")
})

test_that("dgp, allowing for NAs in x and all parameters (loc, scale, shape) w/ log = TRUE", {
  ## Baseline, no missing values in x and any of the parameters
  expect_silent(baseline <- dgp(x = 0:4, loc = 0:4, scale = 1:5, shape = 0:4, log = TRUE))

  ## Missing values in 'x'
  expect_silent(tmp <- dgp(x = c(NA, 1, NA, NA, 4), loc = 0:4, scale = 1:5, shape = 0:4, log = TRUE))
  expect_identical(tmp, baseline * c(NA, 1, NA, NA, 1),
                   info = "expected different return when 'x' includes missing values")

  ## Missing values in 'loc'
  expect_silent(tmp <- dgp(x = 0:4, loc = c(0, NA, 2, NA, NA), scale = 1:5, shape = 0:4, log = TRUE))
  expect_identical(tmp, baseline * c(1, NA, 1, NA, NA),
                   info = "expected different return when 'loc' includes missing values")

  ## Missing values in 'scale'
  expect_silent(tmp <- dgp(x = 0:4, loc = 0:4, scale = c(1, 2, NA, NA, 5), shape = 0:4, log = TRUE))
  expect_identical(tmp, baseline * c(1, 1, NA, NA, 1),
                   info = "expected different return when 'scale' includes missing values")

  ## Missing values in 'shape'
  expect_silent(tmp <- dgp(x = 0:4, loc = 0:4, scale = 1:5, shape = c(NA, 1, 2, 3, NA), log = TRUE))
  expect_identical(tmp, baseline * c(NA, 1, 1, 1, NA),
                   info = "expected different return when 'shape' includes missing values")
})
