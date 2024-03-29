#context("Informative priors: revdbayes vs. evdbayes")

# We check that the the revdbayes (R and C++) functions to evaluate
# informative log-priors agree with each other.

# Set a tolerance for the comparison
my_tol <- 1e-5

# prior_prob

prior_prob_test <- function(x, quant, alpha) {
  r_val <- as.numeric(gev_prob(pars = x, quant = quant, alpha = alpha))
  c_val <- cpp_gev_prob(x = x, ppars = list(quant = quant, alpha = alpha))
  return(c(r_val, c_val))
}

# prior_quant

prior_quant_test <- function(x, prob, shape, scale) {
  r_val <- as.numeric(gev_quant(pars = x, prob = prob, shape = shape,
                                scale = scale))
  c_val <- cpp_gev_quant(x = x, ppars = list(prob = prob, shape = shape,
                                             scale = scale))
  return(c(r_val, c_val))
}

test_function <- function(x, test_string) {
  testthat::test_that(test_string, {
    testthat::expect_equal(x[1], x[2], tolerance = my_tol)
  })
}

test_string <- "prior_prob"

x <- c(84,4.2,-0.3)
quant <- c(85, 88, 95)
alpha <- c(4, 2.5, 2.25, 0.25)
y <- prior_prob_test(x, quant = quant, alpha = alpha)
test_function(y, test_string)

x <- c(84, 1, 0)
y <- prior_prob_test(x, quant = quant, alpha = alpha)
test_function(y, test_string)

test_string <- "prior_quant"

x <- c(43.2, 7.64, 0.32)
prob <- 10^-(1:3)
shape <- c(38.9, 7.1, 47)
scale <- c(1.5, 6.3, 2.6)
y <- prior_quant_test(x, prob = prob, shape = shape, scale = scale)
test_function(y, test_string)

x <- c(50.8, 1.18, 0.65)
y <- prior_quant_test(x, prob = prob, shape = shape, scale = scale)
test_function(y, test_string)
