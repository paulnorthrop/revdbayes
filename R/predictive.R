# ============================ Predictive functions ============================

#' Predictive inference for the largest value observed in N years.
#'
#' Density function, distribution function, quantile function and
#' random generation for the largest value observed in N years using
#' Bayesian (posterior) predictive inference.
#'
#' @param ev_obj An object of class "evpost", a result of a call to
#'   \code{\link{rpost}} with \code{model = "gev"}, \code{model = "os"},
#'   \code{model = "pp"} or \code{model == "bingp"}.  Calling these functions
#'   after a call to \code{rpost} with \code{model == "gp"} will produce an
#'   error, because inferences about the probability of threshold exceedance
#'   are required, in addition to the distribution of threshold excesses.
#'   The model is stored in \code{ev_obj$model}.
#' @param n_years A numeric vector. Values of N.
#' @param npy A numeric scalar. The mean number of observations per year
#'   of data, after excluding any missing values, i.e. the number of
#'   non-missing observations divided by total number of years of non-missing
#'   data.
#'
#' If \code{rpost} was called with \code{model == "bingp"} then \code{npy}
#' must either have been supplied in that call or be supplied here.
#' If \code{npy} is supplied twice then the value supplied here will be
#' used and a warning given.
#'
#' Otherwise, a default value will be assumed if \code{npy} is not supplied,
#' based on the value of \code{model} in the call to \code{rpost}:
#' \itemize{
#'   \item{\code{model = "gev"}:} \code{npy} = 1, i.e. the data were
#'     annual maxima so the block size is one year.
#'   \item{\code{model = "os"}:} \code{npy} = 1, i.e. the data were
#'     annual order statistics so the block size is one year.
#'   \item{\code{model = "pp"}:}
#'     \code{npy} = \code{length(x$data)} / \code{noy},
#'     i.e. the value of \code{noy} used in the call to \code{\link{rpost}}
#'     is equated to a block size of one year.
#' }
#' @param x,q Numeric vectors of quantiles.  If \code{ev_obj$model == "bingp"}
#'   then no element of \code{x} and \code{q} can be less than the
#'   threshold \code{ev_obj$thresh}.
#' @param p A numeric vector of probabilities in (0,1).  If
#'   \code{ev_obj$model == "bingp"} then no element of \code{p} can be less
#'   than the value of \code{ppred(ev_obj, q = ev_obj$thresh)}, i.e.
#'   \code{p} cannot correspond to a predictive quantile that is below the
#'   threshold.
#' @param log A logical scalar.  If TRUE the log-density is returned.
#' @param lower_tail A logical scalar.  If TRUE (default), probabilities
#'   are P[X <= x], otherwise, P[X > x].
#' @details Inferences about future extreme observations are integrated over
#'   the posterior distribution of the model parameters, thereby accounting
#'   for uncertainty in model parameters and uncertainty owing to the
#'   variability of future observations.  In practice the integrals involved
#'   are estimated using an empirical mean over the posterior sample.
#'   See, for example,
#'   \href{http://dx.doi.org/10.1007/978-1-4471-3675-0_9}{Coles (2001),
#'   chapter 9},
#'   \href{http://dx.doi.org/10.1201/b19721-14}{Stephenson (2016)}
#'   or
#'   \href{http://dx.doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   for details.
#'
#'   \strong{GEV / OS / PP}.
#'   If \code{model = "gev"}, \code{model = "os"} or \code{model = "pp"}
#'   in the call to \code{\link{rpost}} then we calculate the number
#'   of blocks in \code{n_years} years, and then convert the posterior
#'   simulated GEV parameters to the \code{n_years} level of aggregation.
#'
#'   In \code{dpred} we calculate using \code{\link{dgev}}
#'   the GEV density at \code{x} for each of the posterior samples
#'   of the location, scale and shape parameters.  Then we take the
#'   mean of these values.
#'
#'   In \code{ppred} we calculate using \code{\link{pgev}}
#'   the GEV density at \code{q} for each of the posterior samples
#'   of the location, scale and shape parameters.  Then we take the
#'   mean of these values.
#'
#'   In \code{qpred} we first calculate initial estimates of the
#'   required predictive quantiles by calculating using \code{\link{qgev}}
#'   the GEV density at \code{p} for each of the posterior samples
#'   of the location, scale and shape parameters and taking the
#'   mean of these values.  Then we solve \code{ppred}(q) = \code{p[i]}
#'   numerically for q for each element \code{p[i]} of \code{p}.
#'
#'   In \code{rpred} for each simulated value of the GEV parameters
#'   at the \code{n_years} level of aggregation we simulate one value from
#'   this GEV distribution using \code{\link{rgev}}.  Thus, each sample
#'   from the predictive distribution is of a size equal to the size of
#'   the posterior sample.
#'
#'   \strong{Binomial-GP}.  If \code{model = "bingp"} in the call to
#'   \code{\link{rpost}} then we calculate the mean number of observations
#'   in \code{n_years} years, i.e. \code{npy * n_years}.
#'
#'   Let \eqn{M_N} be the largest value observed in \eqn{N} years,
#'   \eqn{m} = \code{npy * n_years} and \eqn{u} the threshold
#'   \code{ev_obj$thresh} used in the call to \code{rpost}.
#'   For fixed values of \eqn{\theta = (p, \sigma, \xi)} the distribution
#'   function of \eqn{M_N} is given by \eqn{F(z, \theta)^m}, for
#'   \eqn{z >= u}, where
#'   \deqn{F(z, \theta) = 1 - p * [1 + \xi (x - u) / \sigma] ^ (-1/\xi).}
#'   The distribution function of \eqn{M_N} cannot be evaluated for
#'   \eqn{z < u} because no model has been supposed for observations below
#'   the threshold.
#'
#'   In \code{ppred} we calculate \eqn{F(z, \theta)^m} at \code{q} for
#'   each of the posterior samples \eqn{\theta}.  Then we take the
#'   mean of these values.
#'
#'   In \code{dpred} we calculate the density of of \eqn{M_n}, i.e. the
#'   derivative of \eqn{F(z, \theta)^m} with respect to \eqn{z} at \code{x}
#'   for each of the posterior samples \eqn{\theta}.  Then we take the
#'   mean of these values.
#'
#'   In \code{qpred} we perform a calculation that is analogous to the
#'   GEV case above, i.e. we solve \code{ppred}(q) = \code{p[i]}
#'   numerically for q for each element \code{p[i]} of \code{p}.
#'
#'   In \code{rpred} for each simulated value of the bin-GP parameters
#'   we simulate from the distribution of \eqn{M_N} using the inversion
#'   method applied to the distribution function of \eqn{M_N} given above.
#'   Occasionally a value below the threshold would need to be simulated.
#'   If these instances a missing value code \code{NA} is returned.
#'   Thus, each sample from the predictive distribution is of a size equal
#'   to the size of the posterior sample.
#' @return
#' \itemize{
#'   \item{\code{pred_dgev:}} A \code{length(x)} by \code{length(n_years)}
#'   matrix.  Column i contains the predictive density function
#'   of the \code{n_years[i]} maxiumum evaluated at the values in \code{x}.
#'   \item{\code{pred_pgev:}} A \code{length(q)} by \code{length(n_years)}
#'   matrix.  Column i contains the predictive distribution function
#'   of the \code{n_years[i]} maxiumum evaluated at the values in \code{q}.
#'   \item{\code{pred_qgev:}} A \code{length(p)} by \code{length(n_years)}
#'   matrix.  Column i contains the predictive quantiles of the
#'   \code{n_years[i]} maxiumum evaluated at the probabilities in \code{p} .
#'   \item{\code{pred_rgev:}} An \code{nrow(ev_obj$sim_vals)} by
#'   \code{length(n_years)} matrix.  Column i contains
#'   \code{nrow(ev_obj$sim_vals)} vaues simulated from the predictive
#'   distribution of the \code{n_years[i]} maxiumum.
#' }
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   Chapter 9: \url{http://dx.doi.org/10.1007/978-1-4471-3675-0_9}
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \emph{66}(1), 93-120.
#'   \url{http://dx.doi.org/10.1111/rssc.12159}
#' @references Stephenson, A. (2016). Bayesian Inference for Extreme Value
#'   Modelling. In \emph{Extreme Value Modeling and Risk Analysis: Methods and
#'   Applications}, edited by D. K. Dey and J. Yan, 257-80. London:
#'   Chapman and Hall. \url{http://dx.doi.org/10.1201/b19721-14}
#'   value posterior using the evdbayes package.
#' @examples
#' # GEV
#' data(portpirie)
#' mat <- diag(c(10000, 10000, 100))
#' pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
#' gevp  <- rpost(n = 1000, model = "gev", prior = pn, data = portpirie)
#'
#' q <- seq(4, 7, 0.1)
#' dpred(gevp, x = q)
#' dpred(gevp, q, n_years = c(100, 1000))
#' ppred(gevp, q = q)
#' ppred(gevp, q, n_years = c(100, 1000))
#' p <- c(0.025, 0.25, 0.5, 0.75, 0.975)
#' qq <- qpred(gevp, p)
#' ppred(gevp, qq)
#' qpred(gevp, p, n_years = c(100, 1000))
#'
#' sim1 <- rpred(gevp)
#' sim2 <- pred_rgev(gevp, n_years = c(100, 1000))
#'
#' # Binomial-GP
#' data(gom)
#' u <- quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' bp <- set_bin_prior(prior = "jeffreys")
#' npy_gom <- length(gom)/105
#' gpg1 <- rpost(n = 1000, model = "bingp", prior = fp, thresh = u, data = gom,
#'              bin_prior = bp)
#'
#' q <- seq(10, 30, 0.1)
#' # Setting npy in call to ppred(), for example.
#' ppred(gpg1, q = q, npy = npy_gom)
#'
#' # Setting npy in call to rpost()
#' gpg2 <- rpost(n = 1000, model = "bingp", prior = fp, thresh = u, data = gom,
#'              bin_prior = bp, npy = npy_gom)
#' ppred(gpg2, q = q)
#' dpred(gpg2, x = q)
#' p <- c(0.025, 0.25, 0.5, 0.75, 0.975)
#' qpred(gpg2, p = p)
#'
#' sim1 <- rpred(gpg2)
#' sim2 <- rpred(gpg2, n_years = c(100, 1000))
#' @name predictive_inference
NULL
## NULL

# ----------------------------- dpred ---------------------------------

#' @rdname predictive_inference
#' @export
dpred <- function(ev_obj, x, n_years = 100, npy = NULL, log = FALSE) {
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object of class evpost produced by rpost()")
  }
  if (ev_obj$model == "gp") {
    stop("The model in the call to rpost() cannot be gp.  Use bingp instead.")
  }
  if (ev_obj$model == "bingp") {
    if (is.null(ev_obj$npy) & is.null(npy)) {
      stop("For model=bingp npy must be supplied, here or in call to rpost.")
    }
    if (!is.null(ev_obj$npy) & !is.null(npy)) {
      warning(paste("Two values of npy supplied.  The value npy = ", npy,
                    " from the current call has been used.", sep=""))
    }
    if (!is.null(ev_obj$npy) & is.null(npy)) {
      npy <- ev_obj$npy
    }
  }
  if (ev_obj$model %in% c("gev", "os", "pp")) {
     ret_obj <- pred_dgev(ev_obj = ev_obj, x = x, n_years = n_years,
                          npy = npy, log = log)
  } else if (ev_obj$model == "bingp") {
    ret_obj <- pred_dbingp(ev_obj = ev_obj, x = x, n_years = n_years,
                           npy = npy, log = log)
  } else {
    warning(paste("Predictive functions are not available for model = ",
                  ev_obj$model, sep=""))
    ret_obj <- NULL
  }
  return(ret_obj)
}

# ----------------------------- ppred ---------------------------------

#' @rdname predictive_inference
#' @export
ppred <- function(ev_obj, q, n_years = 100, npy = NULL, lower_tail = TRUE) {
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object of class evpost produced by rpost()")
  }
  if (ev_obj$model == "gp") {
    stop("The model in the call to rpost() cannot be gp.  Use bingp instead.")
  }
  if (ev_obj$model == "bingp") {
    if (is.null(ev_obj$npy) & is.null(npy)) {
      stop("For model=bingp npy must be supplied, here or in call to rpost.")
    }
    if (!is.null(ev_obj$npy) & !is.null(npy)) {
      warning(paste("Two values of npy supplied.  The value npy = ", npy,
                    " from the current call has been used.", sep=""))
    }
    if (!is.null(ev_obj$npy) & is.null(npy)) {
      npy <- ev_obj$npy
    }
  }
  if (ev_obj$model %in% c("gev", "os", "pp")) {
    ret_obj <- pred_pgev(ev_obj = ev_obj, q = q, n_years = n_years,
                         npy = npy, lower_tail = lower_tail)
  } else if (ev_obj$model == "bingp") {
    ret_obj <- pred_pbingp(ev_obj = ev_obj, q = q, n_years = n_years,
                           npy = npy, lower_tail = lower_tail)
  } else {
    warning(paste("Predictive functions are not available for model = ",
                  ev_obj$model, sep=""))
    ret_obj <- NULL
  }
  return(ret_obj)
}

# ----------------------------- qpred ---------------------------------

#' @rdname predictive_inference
#' @export
qpred <- function(ev_obj, p, n_years = 100, npy = NULL, lower_tail = TRUE) {
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object of class evpost produced by rpost()")
  }
  if (ev_obj$model == "gp") {
    stop("The model in the call to rpost() cannot be gp.  Use bingp instead.")
  }
  if (ev_obj$model == "bingp") {
    if (is.null(ev_obj$npy) & is.null(npy)) {
      stop("For model=bingp npy must be supplied, here or in call to rpost.")
    }
    if (!is.null(ev_obj$npy) & !is.null(npy)) {
      warning(paste("Two values of npy supplied.  The value npy = ", npy,
                    " from the current call has been used.", sep=""))
    }
    if (!is.null(ev_obj$npy) & is.null(npy)) {
      npy <- ev_obj$npy
    }
  }
  if (ev_obj$model %in% c("gev", "os", "pp")) {
    ret_obj <- pred_qgev(ev_obj = ev_obj, p = p, n_years = n_years,
                         npy = npy, lower_tail = lower_tail)
  } else if (ev_obj$model == "bingp") {
    ret_obj <- pred_qbingp(ev_obj = ev_obj, p = p, n_years = n_years,
                           npy = npy, lower_tail = lower_tail)
  } else {
    warning(paste("Predictive functions are not available for model = ",
                  ev_obj$model, sep=""))
    ret_obj <- NULL
  }
  return(ret_obj)
}

#' @rdname predictive_inference
#' @export
rpred <- function(ev_obj, n_years = 100, npy = NULL) {
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object of class evpost produced by rpost()")
  }
  if (ev_obj$model == "gp") {
    stop("The model in the call to rpost() cannot be gp.  Use bingp instead.")
  }
  if (ev_obj$model == "bingp") {
    if (is.null(ev_obj$npy) & is.null(npy)) {
      stop("For model=bingp npy must be supplied, here or in call to rpost.")
    }
    if (!is.null(ev_obj$npy) & !is.null(npy)) {
      warning(paste("Two values of npy supplied.  The value npy = ", npy,
                    " from the current call has been used.", sep=""))
    }
    if (!is.null(ev_obj$npy) & is.null(npy)) {
      npy <- ev_obj$npy
    }
  }
  if (ev_obj$model %in% c("gev", "os", "pp")) {
    ret_obj <- pred_rgev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  } else if (ev_obj$model == "bingp") {
    ret_obj <- pred_rbingp(ev_obj = ev_obj, n_years = n_years, npy = npy)
  } else {
    warning(paste("Predictive functions are not available for model = ",
                  ev_obj$model, sep=""))
    ret_obj <- NULL
  }
  return(ret_obj)
}

# ============================ GEV functions ============================

# ----------------------------- pred_dgev ---------------------------------

pred_dgev <- function(ev_obj, x, n_years = 100, npy = NULL, log = FALSE) {
  # Determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  mult <- setup_pred_gev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  #
  loc <- ev_obj$sim_vals[, 1]
  scale <- ev_obj$sim_vals[, 2]
  shape <- ev_obj$sim_vals[, 3]
  n_x <- length(x)
  n_y <- length(n_years)
  d_mat <- matrix(NA, nrow = n_x, ncol = n_y)
  temp <- function(x, loc_n, scale_n, shape) {
    return(mean(dgev(x, loc = loc_n, scale = scale_n, shape = shape)))
  }
  for (i in 1:n_y) {
    # Convert the GEV location and scale parameters from the input time period
    # (block size) to a period of n_years.
    loc_n <- loc + scale * box_cox_vec(x = mult[i], lambda = shape)
    scale_n <- scale * mult[i] ^ shape
    # Calculate the GEV pdf at x for each combination of (loc, scale, shape)
    # in the posterior sample, and take the mean.
    d_mat[, i] <- sapply(x, temp, loc_n = loc_n, scale_n = scale_n,
                         shape = shape)
  }
  if (log) {
    d_mat <- log(d_mat)
  }
  return(d_mat)
}

# ----------------------------- pred_pgev ---------------------------------

pred_pgev <- function(ev_obj, q, n_years = 100, npy = NULL,
                      lower_tail = TRUE) {
  # Determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  mult <- setup_pred_gev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  #
  loc <- ev_obj$sim_vals[, 1]
  scale <- ev_obj$sim_vals[, 2]
  shape <- ev_obj$sim_vals[, 3]
  n_q <- length(q)
  n_y <- length(n_years)
  p_mat <- matrix(NA, nrow = n_q, ncol = n_y)
  temp <- function(q, loc_n, scale_n, shape) {
    return(mean(pgev(q, loc = loc_n, scale = scale_n, shape = shape)))
  }
  for (i in 1:n_y) {
    # Convert the GEV location and scale parameters from the input time period
    # (block size) to a period of n_years.
    loc_n <- loc + scale * box_cox_vec(x = mult[i], lambda = shape)
    scale_n <- scale * mult[i] ^ shape
    # Calculate the GEV cdf at q for each combination of (loc, scale, shape)
    # in the posterior sample, and take the mean.
    p_mat[, i] <- sapply(q, temp, loc_n = loc_n, scale_n = scale_n,
                         shape = shape)
  }
  if (!lower_tail) {
    p_mat <- 1 - p_mat
  }
  return(p_mat)
}

# ----------------------------- pred_qgev ---------------------------------

pred_qgev <- function(ev_obj, p, n_years = 100, npy = NULL,
                      lower_tail = TRUE) {
  # Determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  mult <- setup_pred_gev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  #
  if (!lower_tail) {
    p <- 1 - p
  }
  loc <- ev_obj$sim_vals[, 1]
  scale <- ev_obj$sim_vals[, 2]
  shape <- ev_obj$sim_vals[, 3]
  n_p <- length(p)
  n_y <- length(n_years)
  q_mat <- matrix(NA, nrow = n_p, ncol = n_y)
  temp <- function(p, loc_n, scale_n, shape) {
    return(mean(qgev(p, loc = loc_n, scale = scale_n, shape = shape)))
  }
  for (i in 1:n_y) {
    # Convert the GEV location and scale parameters from the input time period
    # (block size) to a period of n_years.
    loc_n <- loc + scale * box_cox_vec(x = mult[i], lambda = shape)
    scale_n <- scale * mult[i] ^ shape
    # Calculate the GEV quantile at p for each combination of (loc, scale, shape)
    # in the posterior sample, and take the mean.
    #
    # This gives reasonable initial estimates for the predictive quantiles.
    q_mat[, i] <- sapply(p, temp, loc_n = loc_n, scale_n = scale_n,
                         shape = shape)
    #
    logit <- function(p) log(p / (1 - p))
    ob_fn <- function(q, ev_obj, p, n_years, npy) {
      p_val <- pred_pgev(ev_obj = ev_obj, q = q, n_years = n_years, npy = npy)
      (logit(p_val) - logit(p)) ^ 2
    }
    for (j in 1:n_p) {
      qtemp <- stats::nlminb(q_mat[j, i], ob_fn, ev_obj = ev_obj, p = p[j],
                             n_years = n_years[i], npy = npy)
      q_mat[j, i] <- qtemp$par
    }
  }
  return(q_mat)
}

# ----------------------------- pred_rgev ---------------------------------

pred_rgev <- function(ev_obj, n_years = 100, npy = NULL) {
  # Determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  mult <- setup_pred_gev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  #
  loc <- ev_obj$sim_vals[, 1]
  scale <- ev_obj$sim_vals[, 2]
  shape <- ev_obj$sim_vals[, 3]
  n_sim <- length(loc)
  n_y <- length(n_years)
  r_mat <- matrix(NA, nrow = n_sim, ncol = n_y)
  for (i in 1:n_y) {
    # Convert the GEV location and scale parameters from the input time period
    # (block size) to a period of n_years.
    loc_n <- loc + scale * box_cox_vec(x = mult[i], lambda = shape)
    scale_n <- scale * mult[i] ^ shape
    # Simulate a single observation from a GEV distribution corresponding
    # to each parameter combination in the posterior sample.
    r_mat[, i] <- rgev_vec(n = 1, loc = loc_n, scale = scale_n, shape = shape)
  }
  return(r_mat)
}

# -------------------------- setup_pred_gev ------------------------------

setup_pred_gev <- function(ev_obj, n_years, npy) {
  #
  # Determines the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  #
  # Args:
  #   ev_obj  : Object of class evpost return by rpost().
  #   n_years : A numeric vector. Values of N.
  #   npy     : The mean number of observations per year of data, after
  #             excluding any missing values
  #
  # Returns: A numeric scalar.  The value of mult.
  #
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object of class evpost produced by rpost()")
  }
  if (!(ev_obj$model %in% c("gev", "os", "pp"))) {
    stop("The model in the call to rpost() must be gev, os or pp.")
  }
  # If npy is not supplied and the model is GEV or OS then assume that
  # npy = 1, that is, the data were annual maxima or annual order statistics
  # respectively.
  #
  # Also, determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  if (ev_obj$model %in% c("gev", "os")) {
    if (is.null(npy)) {
      npy <- 1
    }
    mult <- n_years * npy
  }
  # Similarly if npy is not supplied and the model is PP then assume that
  # blocks of length one year were set by noy in the call to rpost().
  #
  # Also, determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  if (ev_obj$model == "pp") {
    n <- length(ev_obj$data)
    npy <- ev_obj$noy
    if (is.null(npy)) {
      npy <- n / noy
    }
    mult <- n_years * npy * noy / n
  }
  return(mult)
}

# ============================ binGP functions ============================

# ----------------------------- pred_dbingp ---------------------------------

pred_dbingp <- function(ev_obj, x, n_years = 100, npy = NULL,
                        log = FALSE) {
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object of class evpost produced by rpost()")
  }
  if (ev_obj$model != "bingp") {
    stop("The model in the call to rpost() must be bingp.")
  }
  if (is.null(npy)) {
    stop("npy must be supplied.")
  }
  # Check that q is not less than the threshold used in the call to rpost().
  thresh <- ev_obj$thresh
  if (any(x < thresh)) {
    stop("Invalid x: no element of x can be less than the threshold.")
  }
  # Extract posterior sample of parameters p_u, sigma_u, xi.
  p_u <- ev_obj$bin_sim_vals
  scale <- ev_obj$sim_vals[, 1]
  shape <- ev_obj$sim_vals[, 2]
  # Set up a matrix to contain the results.
  n_x <- length(x)
  n_y <- length(n_years)
  d_mat <- matrix(NA, nrow = n_x, ncol = n_y)
  temp <- function(x, p_u, scale, shape, thresh, mult) {
    # Calculate the distribution function of raw observations, evaluated at q.
    raw_df <- pbingp(q = x, p_u = p_u, loc = thresh, scale = scale,
                     shape = shape)
    # Evaluate the derivative of raw_df ^ mult with respect to x.
    t1 <- mult * exp((mult - 1) * log(raw_df))
    t2 <- p_u * dbingp(x = x, p_u = p_u, loc = thresh, scale = scale,
                       shape = shape)
    # Return the mean of the posterior sample.
    return(mean(t1 * t2))
  }
  # For each value in n_years calculate the distribution function of the
  # n_years maximum.
  mult <- npy * n_years
  for (i in 1:n_y) {
    d_mat[, i] <- sapply(x, temp, p_u = p_u, scale = scale, shape = shape,
                         thresh = thresh, mult = mult[i])
  }
  if (log) {
    d_mat <- log(d_mat)
  }
  return(d_mat)
}

# ----------------------------- pred_pbingp ---------------------------------

pred_pbingp <- function(ev_obj, q, n_years = 100, npy = NULL,
                      lower_tail = TRUE) {
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object of class evpost produced by rpost()")
  }
  if (ev_obj$model != "bingp") {
    stop("The model in the call to rpost() must be bingp.")
  }
  if (is.null(npy)) {
    stop("npy must be supplied.")
  }
  # Check that q is not less than the threshold used in the call to rpost().
  thresh <- ev_obj$thresh
  if (any(q < thresh)) {
    stop("Invalid q: no element of q can be less than the threshold.")
  }
  # Extract posterior sample of parameters p_u, sigma_u, xi.
  p_u <- ev_obj$bin_sim_vals
  scale <- ev_obj$sim_vals[, 1]
  shape <- ev_obj$sim_vals[, 2]
  # Set up a matrix to contain the results.
  n_q <- length(q)
  n_y <- length(n_years)
  p_mat <- matrix(NA, nrow = n_q, ncol = n_y)
  temp <- function(q, p_u, scale, shape, thresh, mult) {
    # Calculate the distribution function of raw observations, evaluated at q.
    raw_df <- pbingp(q = q, p_u = p_u, loc = thresh, scale = scale,
                      shape = shape)
    # Raise this to the power of mult to find the distribution function of
    # the n_year maximum.
    return(mean(exp(mult * log(raw_df))))
  }
  # For each value in n_years calculate the distribution function of the
  # n_years maximum.
  mult <- npy * n_years
  for (i in 1:n_y) {
    p_mat[, i] <- sapply(q, temp, p_u = p_u, scale = scale, shape = shape,
                         thresh = thresh, mult = mult[i])
  }
  if (!lower_tail) {
    p_mat <- 1 - p_mat
  }
  return(p_mat)
}

# ----------------------------- pred_qbingp ---------------------------------

pred_qbingp <- function(ev_obj, p, n_years = 100, npy = NULL,
                        lower_tail = TRUE) {
  if (!lower_tail) {
    p <- 1 - p
  }
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object of class evpost produced by rpost()")
  }
  if (ev_obj$model != "bingp") {
    stop("The model in the call to rpost() must be bingp.")
  }
  if (is.null(npy)) {
    stop("npy must be supplied.")
  }
  # Extract posterior sample of parameters p_u, sigma_u, xi.
  p_u <- ev_obj$bin_sim_vals
  scale <- ev_obj$sim_vals[, 1]
  shape <- ev_obj$sim_vals[, 2]
  # Check that p is not less than the binGP predictive distribution function
  # evaluated at the threshold used in the call to rpost().
  n_y <- length(n_years)
  thresh <- ev_obj$thresh
  p_ok <- rep(TRUE, n_y)
  for (i in 1:n_y) {
    p_min <- pred_pbingp(ev_obj = ev_obj, q = thresh, n_years = n_years[i],
                         npy = npy, lower_tail = lower_tail)
    p_min <- as.vector(p_min)
    if (any(p < p_min)) {
      p_ok[i] <- FALSE
    }
  }
  if (any(!p_ok)) {
    stop("Invalid p for at least one value in n_years")
  }
  n_p <- length(p)
  q_mat <- matrix(NA, nrow = n_p, ncol = n_y)
  temp <- function(p, p_u, loc, scale, shape, mult) {
    pnew <- exp(log(p) / mult)
    return(mean(qbingp(pnew, p_u = p_u, loc = loc, scale = scale,
                       shape = shape)))
  }
  mult <- npy * n_years
  for (i in 1:n_y) {
    # Calculate the quantile at p for each combination of parameters in the
    # in the posterior sample, and take the mean.
    #
    # This gives reasonable initial estimates for the predictive quantiles.
    q_mat[, i] <- sapply(p, temp, p_u = p_u, loc = thresh, scale = scale,
                         shape = shape, mult = mult[i])
    #
    logit <- function(p) log(p / (1 - p))
    ob_fn <- function(q, ev_obj, p, n_years, npy) {
      p_val <- pred_pbingp(ev_obj = ev_obj, q = q, n_years = n_years,
                           npy = npy)
      (logit(p_val) - logit(p)) ^ 2
    }
    for (j in 1:n_p) {
      qtemp <- stats::nlminb(q_mat[j, i], ob_fn, ev_obj = ev_obj, p = p[j],
                             n_years = n_years[i], npy = npy)
      q_mat[j, i] <- qtemp$par
    }
  }
  return(q_mat)
}

pred_rbingp <- function(ev_obj = ev_obj, n_years = n_years, npy = npy) {
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object of class evpost produced by rpost()")
  }
  if (ev_obj$model != "bingp") {
    stop("The model in the call to rpost() must be bingp.")
  }
  if (is.null(npy)) {
    stop("npy must be supplied.")
  }
  # For each value in n_years find the value p_min below which the
  # correspnding simulated value is below the threshold.  A missing value
  # NA is returned to indicate this.
  n_y <- length(n_years)
  thresh <- ev_obj$thresh
  for (i in 1:n_y) {
    p_min <- pred_pbingp(ev_obj = ev_obj, q = thresh, n_years = n_years[i],
                         npy = npy, lower_tail = TRUE)
    p_min <- as.vector(p_min)
  }
  # Extract posterior sample of parameters p_u, sigma_u, xi.
  p_u <- ev_obj$bin_sim_vals
  scale <- ev_obj$sim_vals[, 1]
  shape <- ev_obj$sim_vals[, 2]
  # Set up a matrix to contain the results.
  n_sim <- length(p_u)
  r_mat <- matrix(NA, nrow = n_sim, ncol = n_y)
  #
  mult <- npy * n_years
  # We use the sample underlying simulated uniforms for each value in n_years.
  u <- stats::runif(n_sim)
  for (i in 1:n_y) {
    x_val <- (1 - u ^ (1 / mult[i])) / p_u
    new_mult <- box_cox_vec(x = x_val, lambda = -shape, poly_order = 1)
    r_mat[, i] <- ifelse(u < p_min, NA, thresh - scale * new_mult)
  }
  return(r_mat)
}
