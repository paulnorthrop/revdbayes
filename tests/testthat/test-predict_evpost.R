context("predict.evpost function")

# We check the functions pgev, qgev and rgev.

# Set a tolerance for the comparison of the simulated values
my_tol <- 1e-5

# 1. Check that calling qgev with probabilities p and then calling pgev with
#    the results gets us back to the initial probabilities.

mat <- diag(c(10000, 10000, 100))
pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
gevp  <- rpost_rcpp(n = 1000, model = "gev", prior = pn, data = portpirie)

qs <- predict(gevp, type = "q", n_years = c(100, 1000))$y
ps <- predict(gevp, type = "p", x = qs, n_years = c(100, 1000))$y

check_ps <- matrix(c(0.025, 0.25, 0.5, 0.75, 0.975), 5, 2)

testthat::expect_equal(ps, check_ps, tolerance = my_tol)

