# Add lower bound lb (cf. rprior)
#GEV.prior <- function(GEV.pars,q,N.in,N.out=N.in,log=TRUE,alpha,xi.tol=1e-6,j.max=4){

#' Informative prior for GEV parameters (\eqn{\mu, \sigma, \xi}) constructed on the
#' probability scale.
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param quant A numeric vector of length 3.
#' @param alpha A numeric vector of length 4.
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatability with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gev_prob <- function(pars, quant, alpha, min_xi = -Inf, max_xi = Inf,
                     trendsd = 0){
  # Extract GEV parameter values.
  mu <- pars[1]
  sigma <- pars[2]
  xi <- pars[3]
  # prior is zero if scale parameter is non-positive
  if (sigma <= 0) {
    return(-Inf)
  }
  y <- (quant - mu) / sigma
  lin <- 1 + xi * y
  # prior is zero if any component of lin is not positive.
  if (any(lin <= 0)) {
    return(-Inf)
  }
  #
  # If abs(xi) < xi_tol then xi is treated as being close to zero.
  xi_tol <- 1.e-6
  if (abs(xi) > xi_tol) {
    h <- y / xi - lin * log(lin) / xi ^ 2
  } else {
    j <- 0:4
    h_fun <- function(x) {
      sum((-1) ^ (j+1) * x ^ (j + 2) * xi ^ j / (j + 1) / (j + 2))
    }
    h <- sapply(y, h_fun)
  }
  #
  mat <- matrix(1, 3, 3)
  mat[2, ] <- y
  mat[3, ] <- h
  log_det_mat <- determinant(mat, logarithm = TRUE)$modulus
  # log-density of GEV at q
  log_g <- dgev(x = quant, loc = mu, scale = sigma, shape = xi, log = TRUE)
  pq <- pgev(q = quant, loc = mu, scale = sigma, shape = xi)
  # log-Jacobian
  log_J <- log(sigma) + sum(log_g) + log_det_mat
  # Calculate the Dirichlet log-prior
  # add the endpoints of the probability scale
  pq <- c(0, pq, 1)
  # differences between the pq values
  diff_pq <- diff(pq)
  if (any(diff_pq <= 0)) {
    return(-Inf)
  }
  log_dir_prior <- sum((alpha - 1) * log(diff_pq))
  # Combine the Jacobian with the Dirichlet prior
  val <- log_J + log_dir_prior
  return(val)
}

#' Informative prior for GEV parameters (\eqn{\mu, \sigma, \xi}) constructed on the
#' probability scale.
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param quant A numeric vector of length 3.
#' @param alpha A numeric vector of length 4.
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatability with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gev_quant <- function(pars, prob, shape, scale, min_xi = -Inf, max_xi = Inf,
                      trendsd = 0){
  # Extract GEV parameter values.
  mu <- pars[1]
  sigma <- pars[2]
  xi <- pars[3]
  # prior is zero if scale parameter is non-positive
  if (sigma <= 0) {
    return(-Inf)
  }
  # Calculate quantiles. Note: prob contains exceedance probabilities
  quant <- qgev(p = 1 - prob, loc = mu, scale = sigma, shape = xi)
  y <- (quant - mu) / sigma
  lin <- 1 + xi * y
  # prior is zero if any component of lin is not positive.
  if (any(lin <= 0)) {
    return(-Inf)
  }
  #
  # If abs(xi) < xi_tol then xi is treated as being close to zero.
  xi_tol <- 1.e-6
  if (abs(xi) > xi_tol) {
    h <- y / xi - lin * log(lin) / xi ^ 2
  } else {
    j <- 0:4
    h_fun <- function(x) {
      sum((-1) ^ (j+1) * x ^ (j + 2) * xi ^ j / (j + 1) / (j + 2))
    }
    h <- sapply(y, h_fun)
  }
  #
  mat <- matrix(1, 3, 3)
  mat[2, ] <- y
  mat[3, ] <- h
  log_det_mat <- determinant(mat, logarithm = TRUE)$modulus
  # log-Jacobian
  log_J <- log(sigma) + log_det_mat
  # Calculate the gamma log-prior
  # differences between the quant values
  diff_quant <- diff(c(0, quant))
  log_gamma_prior <- sum((shape - 1) * log(diff_quant) - diff_quant / scale)
  # Combine the Jacobian with the Dirichlet prior
  val <- log_J + log_gamma_prior
  return(val)
}
