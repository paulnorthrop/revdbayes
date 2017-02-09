# ================================ rprior_quant ===============================

#' Prior simulation of GEV parameters - prior on quantile scale
#'
#' Simulates from the prior distribution for GEV parameters proposed in
#' Coles and Tawn (1996), based on independent gamma priors for differences
#' between quantiles.
#'
#' @param n A numeric scalar. The size of posterior sample required.
#' @param prob A numeric vector of length 3. Exceedance probabilities
#'   corresponding to the quantiles used to specify the prior distribution.
#' @param shape A numeric vector of length 3. Respective shape parameters of
#'   the gamma priors for the quantile differences.
#' @param scale A numeric vector of length 3. Respective scale parameters of
#'   the gamma priors for the quantile differences.
#'
#' @details The simulation is based on the way that the prior is constructed.
#'   [See Coles and Tawn (1996), Stephenson (1996) or the evdbayes user guide
#'   for details of the construction of the prior.] First, the quantile
#'   differences are simulated from the specified gamma distributions.
#'   Then the simulated quantiles are calculated. Then the GEV location,
#'   scale and shape parameters that give these quantile values are found,
#'   by solving numerically a set of three non-linear equations in which the
#'   GEV quantile function is equated to the simulated quantiles.  This is
#'   reduced to a one-dimensional optimisation over the GEV shape parameter.
#' @return An \code{n} by 3 numeric matrix.
#' @seealso \code{\link[evdbayes]{prior.quant}} to set this prior using the
#'   evdbayes package.
#' @seealso \code{\link[evdbayes]{posterior}}: evdbayes function that can sample from the
#'   prior distribution (using MCMC) if the argument \code{lh = "none"} is
#'   given.
#' @seealso \code{\link{rpost}} for sampling from an extreme value posterior
#'   distribution.
#' @references Coles, S. G. and Tawn, J. A. (1996) A Bayesian analysis of
#'   extreme rainfall data. \emph{Appl. Statist.}, \strong{45}, 463-478.
#'   \url{http://dx.doi.org/10.2307/2986068}.
#' @references Stephenson, A. 2016. Bayesian Inference for Extreme Value
#'   Modelling. In \emph{Extreme Value Modeling and Risk Analysis: Methods and
#'   Applications}, edited by D. K. Dey and J. Yan, 257-80. London:
#'   Chapman and Hall. \url{http://dx.doi.org/10.1201/b19721-14}
#'
#' @examples
#' my_sh <- c(38.9,7.1,47)
#' my_sc <- c(1.5,6.3,2.6)
#' my_pr <- 10^-(1:3)
#' x1 <- rprior_quant(n = 1000, shape = my_sh, scale = my_sc, prob = my_pr)
#' @export
rprior_quant <- function(n, shape, scale, prob){
  #
  pq <- prob
  #
  # Simulate qp_tildes
  #
  qp1_tilde <- rgamma(n, shape = shape[1], scale = scale[1])
  qp2_tilde <- rgamma(n, shape = shape[2], scale = scale[2])
  qp3_tilde <- rgamma(n, shape = shape[3], scale = scale[3])
  #
  # Transform to qps
  #
  qp1 <- qp1_tilde
  qp2 <- qp1_tilde + qp2_tilde
  qp3 <- qp1_tilde + qp2_tilde + qp3_tilde
  qp.sim <- cbind(qp1, qp2, qp3)
  #
  # Transform to (mu, sigma, xi)
  #
  f.xi <- function(xi, pq){
    xp <- -log(1 - pq)
    ifelse(abs(xi) < 1e-6, -log(xp) + xi * log(xp) ^ 2 / 2,
           (xp ^ (-xi) - 1) / xi)
  }
  #
  obfn <- function(xi, qp, pq){
    f1 <- f.xi(xi = xi, pq = pq[1])
    f2 <- f.xi(xi = xi, pq = pq[2])
    f3 <- f.xi(xi = xi, pq = pq[3])
    sigma <- (qp[3] - qp[1])/(f3 - f1)
    mu <- qp[1] - sigma * f1
    (mu + sigma * f2 - qp[2])^2
  }
  #
  pars <- matrix(NA, ncol=3, nrow=n)
  for (i in 1:n){
    xi <- nlminb(0.1, obfn, qp = qp.sim[i, ], pq = pq)$par
    f1 <- f.xi(xi =xi, pq = pq[1])
    f3 <- f.xi(xi =xi, pq = pq[3])
    sigma <- (qp.sim[i,3] - qp.sim[i, 1])/(f3 - f1)
    mu <- qp.sim[i, 1] - sigma * f1
    pars[i,] <- c(mu, sigma, xi)
  }
  #
  return(pars)
}
