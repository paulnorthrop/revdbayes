# ================================= kgaps_post ================================

#' Random sampling from K-gaps posterior distribution
#'
#' Uses the \code{\link[rust]{rust}} package to simulate from the posterior
#' distribution of the extremal index \eqn{\theta} based on the K-gaps model
#' for threshold interexceedance times of Suveges and Davison (2010).
#'
#' @param n A numeric scalar. The size of posterior sample required.
#' @param data A numeric vector or numeric matrix of raw data.  If \code{data}
#'   is a matrix then the log-likelihood is constructed as the sum of
#'   (independent) contributions from different columns. A common situation is
#'   where each column relates to a different year.
#'
#'   If \code{data} contains missing values then
#'   \code{\link[exdex]{split_by_NAs}} is used to divide the data further into
#'   sequences of non-missing values, stored in different columns in a matrix.
#'   Again, the log-likelihood is constructed as a sum of contributions from
#'   different columns.
#' @param thresh A numeric scalar.  Extreme value threshold applied to data.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = \max(T - K, 0)}{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from right-censored inter-exceedance times, relating to the
#'   first and last observations.  It is known that these times are greater
#'   than or equal to the time observed.
#'   If \code{data} has multiple columns then there will be right-censored
#'   first and last inter-exceedance times for each column.  See also the
#'   \strong{Details} section of \code{\link[exdex]{kgaps}}.
#' @param alpha,beta Positive numeric scalars.  Parameters of a
#'   beta(\eqn{\alpha}, \eqn{\beta}) prior for \eqn{\theta}.
#' @param param A character scalar.  If \code{param = "logit"} (the default)
#'   then we simulate from the posterior distribution of
#'   \eqn{\phi = \log(\theta / (1-\theta))}{\phi = log(\theta / (1-\theta))}
#'   and then transform back to the
#'   \eqn{\theta}-scale.  If \code{param = "theta"} then we simulate
#'   directly from the posterior distribution of \eqn{\theta}, unless
#'   the sample K-gaps are all equal to zero or all positive, when we revert
#'   to \code{param = "logit"}.  This is to avoid sampling directly from a
#'   posterior with mode equal to 0 or 1.
#' @param use_rcpp A logical scalar.  If \code{TRUE} (the default) the
#'   rust function \code{\link[rust]{ru_rcpp}} is used for
#'   posterior simulation.  If \code{FALSE} the (slower) function
#'   \code{\link[rust]{ru}} is used.
#' @details A beta(\eqn{\alpha}, \eqn{\beta}) prior distribution is used for
#'   \eqn{\theta} so that the posterior from which values are simulated is
#'   proportional to
#'   \deqn{\theta ^ {2 N_1 + \alpha - 1} (1 - \theta) ^ {N_0 + \beta - 1}
#'     \exp\{- \theta q (S_0 + \cdots + S_N)\}.}{%
#'     \theta ^ (2 N_1 + \alpha - 1) * (1 - \theta) ^ (N_0 + \beta - 1) *
#'     exp(- \theta q (S_0 + \dots + S_N)).}
#'   See \code{\link[exdex]{kgaps_stat}} for a description of the variables
#'   involved in the contribution of the likelihood to this expression.
#'
#'   The \code{\link[rust]{ru}} function in the \code{\link[rust]{rust}}
#'   package simulates from this posterior distribution using the
#'   generalised ratio-of-uniforms distribution.  To improve the probability
#'   of acceptance, and to ensure that the simulation will work even in
#'   extreme cases where the posterior density of \eqn{\theta} is unbounded as
#'   \eqn{\theta} approaches 0 or 1, we simulate from the posterior
#'   distribution of
#'   \eqn{\phi = \log(\theta / (1-\theta))}{\phi = log(\theta / (1-\theta))}
#'   and then transform back to the \eqn{\theta}-scale.
#' @return An object (list) of class \code{"evpost"}, which has the same
#'   structure as an object of class \code{"ru"} returned from
#'   \code{\link[rust]{ru}}.
#'   In addition this list contains
#'   \itemize{
#'     \item{\code{call}:} The call to \code{kgaps()}.
#'     \item{\code{model}:} The character scalar \code{"kgaps"}.
#'     \item{\code{thresh}:} The argument \code{thresh}.
#'     \item{\code{ss}:} The sufficient statistics for the K-gaps likelihood,
#'       as calculated by \code{\link[exdex]{kgaps_stat}}.
#'   }
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221. \doi{10.1214/09-AOAS292}
#' @seealso \code{\link[rust]{ru}} for the form of the object returned by
#'   \code{kgaps_post}.
#' @seealso \code{\link{dgaps_post}} for Bayesian inference about the
#'   extremal index \eqn{\theta} using the \eqn{D}-gaps model.
#' @examples
#' # Newlyn sea surges
#'
#' thresh <- quantile(newlyn, probs = 0.90)
#' k_postsim <- kgaps_post(newlyn, thresh)
#' plot(k_postsim)
#'
#' ### Cheeseboro wind gusts
#'
#' k_postsim <- kgaps_post(exdex::cheeseboro, thresh = 45, k = 3)
#' plot(k_postsim)
#' @export
kgaps_post <- function(data, thresh, k = 1, n = 1000, inc_cens = TRUE,
                       alpha = 1, beta = 1, param = c("logit", "theta"),
                       use_rcpp = TRUE) {
  if (!is.numeric(thresh) || length(thresh) != 1) {
    stop("thresh must be a numeric scalar")
  }
  if (thresh >= max(data, na.rm = TRUE)) {
    stop("thresh must be less than max(data)")
  }
  if (!is.numeric(k) || k <= 0 || length(k) != 1) {
    stop("k must be a positive scalar")
  }
  if (alpha <= 0 | beta <= 0) {
    stop("alpha and beta must be positive.")
  }
  param <- match.arg(param)
  # Use exdex::kgaps() to find the MLE and get the sufficient statistics
  temp <- exdex::kgaps(data = data, u = thresh, k = k, inc_cens = inc_cens)
  theta_mle <- temp$theta
  ss <- temp$ss
  #
  # Set an initial value for theta, and perhaps phi. We do this by noting that
  # (a) the Jacobian of the transformation from theta to phi multiplies the
  # likelihood by a factor of theta * (1 - theta), and (b) the beta prior
  # includes a factor of theta ^ (alpha - 1) * (1 - theta) ^ (beta - 1).
  # Therefore, we can find the MAP of theta by solving a quadratic equation
  # for theta, in the same way that we do to find the MLE for theta.
  #
  theta_init <- kgaps_quad_solve(ss$N0 + alpha, ss$N1 + beta, ss$sum_qs)
  # Set essential arguments to ru()
  # ... including the K-gaps posterior distribution for theta
  if (use_rcpp) {
    post_ptr <- kgaps_logpost_xptr("kgaps")
    for_post <- c(ss, list(alpha = alpha, beta = beta))
    for_ru <- list(logf = post_ptr, pars = for_post, n = n)
  } else {
    logpost <- function(theta, ss) {
      loglik <- do.call(kgaps_loglik, c(list(theta = theta), ss))
      if (is.infinite(loglik)) return(loglik)
      # Add beta(alpha, beta) prior
      logprior <- stats::dbeta(theta, alpha, beta, log = TRUE)
      return(loglik + logprior)
    }
    for_ru <- list(logf = logpost, ss = ss, n = n)
  }
  # Sample on the logit scale phi = log(theta / (1 - theta)) ?
  if (ss$N0 == 0 || ss$N1 == 0) {
    param <- "logit"
  }
  if (param == "logit") {
    # Transformation, Jacobian and initial estimate
    if (use_rcpp) {
      phi_to_theta <- phi_to_theta_xptr("kgaps")
      log_j <- log_j_xptr("kgaps")
    } else {
      phi_to_theta <- function(phi) {
        ephi <- exp(phi)
        return(ephi / (1 + ephi))
      }
      log_j <- function(theta) {
        return(-log(theta) - log(1 - theta))
      }
    }
    phi_init <- log(theta_init / (1 - theta_init))
    trans_list <- list(phi_to_theta = phi_to_theta, log_j = log_j)
    for_ru <- c(for_ru, list(init = phi_init, trans = "user"), trans_list)
  } else {
    for_ru <- c(for_ru, list(init = theta_init, trans = "none"))
  }
  ru_fn <- ifelse(use_rcpp, rust::ru_rcpp, rust::ru)
  temp <- do.call(ru_fn, for_ru)
  temp$model <- "kgaps"
  temp$thresh <- thresh
  temp$ss <- ss
  colnames(temp$sim_vals) <- "theta"
  temp$call <- match.call(expand.dots = TRUE)
  class(temp) <- "evpost"
  return(temp)
}

# =============================== kgaps_loglik ================================

kgaps_loglik <- function(theta, N0, N1, sum_qs, n_kgaps){
  if (theta < 0 || theta > 1) {
    return(-Inf)
  }
  loglik <- 0
  if (N1 > 0) {
    loglik <- loglik + 2 * N1 * log(theta) - sum_qs * theta
  }
  if (N0 > 0) {
    loglik <- loglik + N0 * log(1 - theta)
  }
  return(loglik)
}

# ============================== kgaps_quad_solve =============================

kgaps_quad_solve <- function(N0, N1, sum_qs) {
  aa <- sum_qs
  bb <- -(N0 + 2 * N1 + sum_qs)
  cc <- 2 * N1
  qq <- -(bb - sqrt(bb ^ 2 - 4 * aa * cc)) / 2
  theta_mle <- cc / qq
  return(theta_mle)
}
