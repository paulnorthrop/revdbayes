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
  # If abs(xi) < x_tol then xi is treated as being close to zero.
  xi_tol <- 1.e-6
  dd <- 2 * xi_tol
  j_max <- 4
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
  if (abs(xi) > xi_tol) {
    h <- y / xi - lin * log(lin) / xi ^ 2
  } else {
    j <- 0:j_max
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
  # Combine the Jacobian with the Dirichlet prior
  val <- log_J + Dir_log_prior(pq, alpha = alpha)
  return(val)
}

Dir_log_prior <- function(pq, alpha = NULL){
  #
  # pq    : 3-vector of reference non-exceedance probabilities (i.e cdf)
  # alpha : (length(pq)+1)-vector of Dirichlet parameters
  #
  if (is.null(alpha)) {
    alpha <- rep(1, length(pq) + 1)
  }
  # add the endpoints of the probability scale
  pq <- c(0, pq, 1)
  # differences between the pq values
  diff_pq <- diff(pq)
  if (any(diff_pq <= 0)) {
    return(-Inf)
  }
  # log-prior, up to an additive constant
  return(sum((alpha - 1) * log(diff_pq)))
}
