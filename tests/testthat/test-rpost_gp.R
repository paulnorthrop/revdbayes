context("GP: rpost_rcpp vs rpost")

# We check that the values simulated using rpost() and rpost_rcpp() are
# (close enough to) identical when they are called using the same data,
# the same prior and starting from the same random number seed.

# Set a tolerance for the comparison of the simulated values
my_tol <- 1e-5

# Generalized Pareto distribution --------------------------

#model, prior, rotate, trans, use_phi_map

# Set up the data and the threshold.
data(gom)
u <- quantile(gom, probs = 0.65)

gp_test <- function(seed = 47, prior, model = "gp", n = 10, rotate = TRUE,
                    trans = "none", use_phi_map = FALSE, ...){
  if (prior == "user") {
    prior_rfn <- set_prior(prior = "flat", model = "gp", min_xi = -1)
    ptr_gp_flat <- create_prior_xptr("gp_flat")
    prior_cfn <- set_prior(prior = ptr_gp_flat, model = "gp", min_xi = -1)
  } else {
    prior_rfn <- set_prior(prior = prior, model = model, ...)
    prior_cfn <- set_prior(prior = prior, model = model, ...)
  }
  set.seed(seed)
  res1 <- rpost(n = n, model = "gp", prior = prior_rfn, thresh = u,
                     data = gom, rotate = rotate, trans = trans,
                     use_phi_map = use_phi_map)
  set.seed(seed)
  res2 <- rpost_rcpp(n = n, model = "gp", prior = prior_cfn, thresh = u,
                     data = gom, rotate = rotate, trans = trans,
                     use_phi_map = use_phi_map)
  return(list(sim1 = res1$sim_vals, sim2 = res2$sim_vals))
}

test_function <- function(x, test_string) {
  testthat::test_that(test_string, {
    skip_on_cran()
    testthat::expect_equal(x$sim1, x$sim2, tolerance = my_tol)
  })
}

rotate_vals <- c(FALSE, TRUE)
trans_vals <- c("none", "BC")
use_phi_map_vals <- c(FALSE, TRUE)

for (rotate in rotate_vals) {
  for (trans in trans_vals) {
    for (use_phi_map in use_phi_map_vals) {
      test_string <- paste("rotate =", rotate, "trans =", trans,
                            "use_phi_map =", use_phi_map)
      x <- gp_test(prior = "flat", min_xi = -1,
                   rotate = rotate, trans = trans, use_phi_map = use_phi_map)
      test_function(x, test_string)
      x <- gp_test(prior = "flatflat", min_xi = -1,
                   rotate = rotate, trans = trans, use_phi_map = use_phi_map)
      test_function(x, test_string)
      x <- gp_test(prior = "norm", mean = c(0,0), cov = diag(100,2),
                   rotate = rotate, trans = trans, use_phi_map = use_phi_map)
      test_function(x, test_string)
      x <- gp_test(prior = "mdi",
                   rotate = rotate, trans = trans, use_phi_map = use_phi_map)
      test_function(x, test_string)
      x <- gp_test(prior = "jeffreys",
                   rotate = rotate, trans = trans, use_phi_map = use_phi_map)
      test_function(x, test_string)
      x <- gp_test(prior = "beta",
                   rotate = rotate, trans = trans, use_phi_map = use_phi_map)
      test_function(x, test_string)
      x <- gp_test(prior = "user",
                   rotate = rotate, trans = trans, use_phi_map = use_phi_map)
      test_function(x, test_string)
    }
  }
}