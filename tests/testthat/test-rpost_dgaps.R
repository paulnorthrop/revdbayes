#context("D-gaps posterior: R vs Rcpp")

# We check that the values simulated using dgaps_post are (close enough to)
# identical when they are called using use_rcpp = TRUE and use_rcpp = FALSE
# using the same prior and starting from the same random number seed.

dgaps_test <- function(seed = 47, data = newlyn, thresh, D = 1, n = 5,
                       inc_cens = TRUE, alpha = 1, beta = 1, param = "logit"){
  set.seed(seed)
  res1 <- dgaps_post(data = data, thresh = thresh, D = D, n = n,
                     inc_cens = inc_cens, alpha = alpha, beta = beta,
                     param = param, use_rcpp = TRUE)
  set.seed(seed)
  res2 <- dgaps_post(data = data, thresh = thresh, D = D, n = n,
                     inc_cens = inc_cens, alpha = alpha, beta = beta,
                     param = param, use_rcpp = FALSE)
  return(list(sim1 = as.numeric(res1$sim_vals),
              sim2 = as.numeric(res2$sim_vals)))
}

# Set a tolerance for the comparison of the simulated values
my_tol <- 1e-5

test_function <- function(x, test_string) {
  testthat::test_that(test_string, {
    testthat::expect_equal(x$sim1, x$sim2, tolerance = my_tol)
  })
}

# Set a threshold
thresh <- stats::quantile(newlyn, probs = 0.90)

x <- dgaps_test(data = newlyn, thresh = thresh)
test_function(x, "inc_cens = TRUE, param = logit")

x <- dgaps_test(data = newlyn, thresh = thresh, inc_cens = FALSE)
test_function(x, "inc_cens = FALSE, param = logit")

x <- dgaps_test(data = newlyn, thresh = thresh, param = "theta")
test_function(x, "inc_cens = TRUE, param = theta")

x <- dgaps_test(data = newlyn, thresh = thresh, inc_cens = FALSE,
                param = "theta")
test_function(x, "inc_cens = FALSE, param = theta")
