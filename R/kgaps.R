# ================================= kgaps_mle =================================
#
#' Maximum likelihood estimation for the K-gaps model
#'
#' Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
#' based on the K-gaps model of Suveges and Davison (2010).
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.
#' @param conf  A numeric scalar.  If \code{conf} is supplied then a
#'   \code{conf}\% likelihood-based confidence interval for \eqn{\theta} is
#'   estimated.
#' @details The maximum likelihood estimate of the extremal index \eqn{\theta}
#'   under the K-gaps model of Suveges and Davison (2010) is calculated.
#'   If \code{inc_cens = TRUE} then information from censored inter-exceedance
#'   times is included in the likelihood to be maximised, following
#'   Attalides (2015).  The form of the log-likelihood is given in the
#'   \strong{Details} section of \code{\link{kgaps_stats}}.
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{http://dx.doi.org/10.1214/09-AOAS292"}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{http://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @return A list containing
#'   \itemize{
#'     \item {\code{theta_mle} : } {The maximum likelihood estimate (MLE) of
#'       \eqn{\theta}.}
#'     \item {\code{theta_se} : } {The estimated standard error of the MLE.}
#'     \item {\code{conf_int} : } {(If \code{conf} is supplied) a numeric
#'       vector of length two giving lower and upper confidence limits for
#'       \eqn{\theta}.}
#'   }
#' @seealso \code{\link{kgaps_stats}} for the calculation of sufficient
#'   statistics for the K-gaps model.
#' @examples
#' u <- quantile(newlyn, probs = 0.90)
#' # MLE and SE only
#' kgaps_mle(newlyn, u)
#' # MLE, SE and 95% confidence interval
#' kgaps_mle(newlyn, u, conf = 95)
kgaps_mle <- function(data, u, k = 1, inc_cens = FALSE, conf = NULL) {
  # Calculate sufficient statistics
  ss <- kgaps_stats(data, u, k, inc_cens)
  # If N0 = 0 then all exceedances occur singly (all K-gaps are positive)
  # and the likelihood is maximised at theta = 1.
  N0 <- ss$N0
  # If N1 = 0 then we are in the degenerate case where there is one cluster
  # (all K-gaps are zero) and the likelihood is maximised at theta = 0.
  N1 <- ss$N1
  if (N1 == 0) {
    theta_mle <- 0
  } else if (N0 == 0) {
    theta_mle <- 1
  } else {
    sum_qs <- ss$sum_qs
    aa <- sum_qs
    bb <- -(N0 + 2 * N1 + sum_qs)
    cc <- 2 * N1
    qq <- -(bb - sqrt(bb ^ 2 - 4 * aa * cc)) / 2
    theta_mle <- cc / qq
  }
  # Estimate standard error
  obs_info <- 0
  if (N0 > 0) {
    obs_info <- obs_info + N0 / (1 - theta_mle) ^ 2
  }
  if (N1 > 0) {
    obs_info <- obs_info + 2 * N1 / theta_mle ^ 2
  }
  theta_se <- sqrt(1 / obs_info)
  if (is.null(conf)) {
    return(list(theta_mle = theta_mle, theta_se = theta_se))
  }
  conf_int <- kgaps_conf_int(theta_mle, ss, conf = 95)
  return(list(theta_mle = theta_mle, theta_se = theta_se, conf_int = conf_int))
}

# ================================ kgaps_stats ================================

#' Sufficient statistics for the K-gaps model
#'
#' Calculates sufficient statistics for the K-gaps model for the extremal index
#' \eqn{\theta}.
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.
#' @details The sample K-gaps are \eqn{S_0, S_1, ..., S_(N-1),
#'   S_N}, where \eqn{S_1, ..., S_(N-1)} are uncensored and \eqn{S_0} and
#'   \eqn{S_N} are censored.  Under the assumption that the K-gaps are
#'   independent, the log-likelihood of the K-gaps model is given by
#'   \deqn{l(\theta; S_0, ..., S_N) = N_0 log(1 - \theta) + 2 N_1 log \theta -
#'     \theta q (S_0 + ... + S_N),}
#'    where \eqn{q} is the threshold exceedance probability,
#'    \eqn{N_0} is the number of sample K-gaps that are equal to zero and
#'    (apart from an adjustment for the contributions of \eqn{S_0} and
#'    \eqn{S_N}) \eqn{N_1} is the number of positive sample K-gaps.
#'    Specifically, \eqn{N_1} is equal to the number of \eqn{S_1, ..., S_(N-1)}
#'    that are positive plus \eqn{(I_0 + I_N) / 2}, where \eqn{I_0 = 1} if
#'    \eqn{S_0} is greater than zero and similarly for \eqn{I_N}.
#'    The differing treatment of uncensored and censored K-gaps reflects
#'    differing contributions to the likelihood.
#'    For full details see Suveges and Davison (2010) and Attalides (2015).
#' @return A list containing the sufficient statistics, with components
#'   \itemize{
#'     \item {\code{N0} : } {the number of zero K-gaps}
#'     \item {\code{N1} : } {contribution from non-zero K-gaps (see
#'       \strong{Details})}
#'     \item {\code{sum_qs} : } {the sum of the (scaled) K-gaps, i.e.
#'       \eqn{q (S_0 + ... + S_N)}, where \eqn{q} is estimated by the
#'       proportion of threshold exceedances.}
#'   }
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{http://dx.doi.org/10.1214/09-AOAS292"}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{http://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @seealso \code{\link{kgaps_mle}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @examples
#' u <- quantile(newlyn, probs = 0.90)
#' kgaps_stats(newlyn, u)
kgaps_stats <- function(data, u, k = 1, inc_cens = FALSE) {
  if (any(is.na(data))) {
    stop("No missing values are allowed in ''data''")
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > u]
  N_u <- length(exc_u)
  q_u <- N_u / nx
  # Inter-exceedances times and K-gaps
  T_u <- diff(exc_u)
  S_k <- pmax(T_u - k,0)
  # N0, N1, sum of scaled K-gaps
  N1 <- sum(S_k > 0)
  N0 <- N_u - 1 - N1
  sum_qs <- sum(q_u * S_k)
  # Include censored inter-exceedance times?
  if (inc_cens) {
    # censored inter-exceedance times and K-gaps
    T_u_cens <- c(exc_u[1] - 1, nx - exc_u[N_u])
    S_k_cens <- pmax(T_u_cens - k, 0)
    # N0, N1, sum of scaled K-gaps
    N1_cens <- sum(S_k_cens > 0)
    sum_s_cens <- sum(q_u * S_k_cens)
    # Add contributions.
    # Note: we divide N1_cens by two because a censored non-zero K-gap S_c
    # contributes theta exp(-theta q_u S_c) to the K-gaps likelihood,
    # whereas a non-censored non-zero K-gap contributes
    # theta^2 exp(-theta q_u S_c).
    # See equation (4.3) of Attalides (2015)
    N1 <- N1 + N1_cens / 2
    sum_qs <- sum_qs + sum_s_cens
  }
  return(list(N0 = N0, N1 = N1, sum_qs = sum_qs))
}

# =============================== kgaps_loglik ================================

kgaps_loglik <- function(theta, N0, N1, sum_qs){
  loglik <- 0
  if (N1 > 0) {
    loglik <- loglik + 2 * N1 * log(theta) - sum_qs * theta
  }
  if (N0 > 0) {
    loglik <- loglik + N0 * log(1 - theta)
  }
  return(loglik)
}

# ============================== kgaps_conf_int ===============================

kgaps_conf_int <- function(theta_mle, ss, conf = 95) {
  cutoff <- qchisq(conf / 100, df = 1)
  theta_list <- c(list(theta = theta_mle), ss)
  max_loglik <- do.call(kgaps_loglik, theta_list)
  ob_fn <- function(theta) {
    theta_list$theta <- theta
    loglik <- do.call(kgaps_loglik, theta_list)
    return(2 * (max_loglik - loglik) - cutoff)
  }
  ci_low <- 0
  ci_up <- 1
  if (ss$N1 > 0) {
    ci_low <- uniroot(ob_fn, c(0, theta_mle))$root
  }
  if (ss$N0 > 0) {
    ci_up <- uniroot(ob_fn, c(theta_mle, 1))$root
  }
  return(c(ci_low, ci_up))
}
