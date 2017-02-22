# ============================ Predictive GEV functions ============================

#' GEV-based predictive inference for the largest value observed in N years.
#'
#' Density function, distribution function, quantile function and
#' random generation for the largest value observed in N years using
#' Bayesian predictive inference based on a generalized extreme value (GEV)
#' distribution.
#'
#' @param ev_obj An object of class "evpost", a result of a call to
#'   \code{\link{rpost}} with \code{model = "gev"}, \code{model = "os"}
#'   or \code{model = "pp"}.
#' @param n_years A numeric vector. Values of N.
#' @param npy The mean number of observations per year of data, after
#'   excluding any missing values, i.e. the number of non-missing observations
#'   divided by total number of years of non-missing data.
#'
#' A default value will be assumed if \code{npy} is not supplied:
#' \itemize{
#'   \item{\code{model = "gev"}:} \code{npy} = 1, i.e. the data were
#'     annual maxima.
#'   \item{\code{model = "os"}:} \code{npy} = 1, i.e. the data were
#'     annual order statistics.
#'   \item{\code{model = "pp"}:}
#'     \code{npy} = \code{length(x$data)} / \code{noy},
#'     i.e. the value of \code{noy} used in the call to \code{\link{rpost}}
#'     equated to a block size of one year.
#' }
#' @param x,q Numeric vectors of quantiles.
#' @param p A numeric vector of probabilities in (0,1).
#' @param log A logical scalar.  If TRUE the log-density is returned.
#' @param lower_tail A logical scalar.  If TRUE (default), probabilities
#'   are P[X <= x], otherwise, P[X > x].
#' @details In each function we first infer the number, mult, of blocks in
#'   n_years years, and then convert the posterior simulated GEV
#'   parameters to the \code{n_years} level of aggregation.
#'
#'   In \code{pred_dgev} we calculate using \code{\link{dgev}}
#'   the GEV density at \code{x} for each of the posterior samples
#'   of the location, scale and shape parameters.  Then we take the
#'   mean of these values.
#'
#'   In \code{pred_pgev} we calculate using \code{\link{pgev}}
#'   the GEV density at \code{q} for each of the posterior samples
#'   of the location, scale and shape parameters.  Then we take the
#'   mean of these values.
#'
#'   In \code{pred_qgev} we first calculate initial estimates of the
#'   required predictive quantiles: we calculate using \code{\link{qgev}}
#'   the GEV density at \code{p} for each of the posterior samples
#'   of the location, scale and shape parameters.  Then we take the
#'   mean of these values.  Then we solve \code{pred_pgev}(q) = \code{p[i]}
#'   numerically for q for each element \code{p[i]} of \code{p}.
#'
#'   In \code{pred_rgev} for each simulated value of the GEV parameters
#'   at the \code{n_years} we simulate one value from this GEV distribution
#'   using \code{\link{rgev}}.  Thus, each sample from the predictive
#'   distribution is of a size equal to the size of the posterior sample.
#'
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
#' @examples
#' data(portpirie)
#' mat <- diag(c(10000, 10000, 100))
#' pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
#' gevp  <- rpost(n = 1000, model = "gev", prior = pn, data = portpirie)
#'
#' pred_dgev(gevp, q)
#' pred_dgev(gevp, q, n_years <- c(100, 1000))
#'
#' q <- seq(4, 7, 0.1)
#' pred_pgev(gevp, q)
#' pred_pgev(gevp, q, n_years <- c(100, 1000))
#'
#' p <- c(0.025, 0.25, 0.5, 0.75, 0.975)
#' qq <- pred_qgev(gevp, p)
#' pred_pgev(gevp, qq)
#' pred_qgev(gevp, p, n_years <- c(100, 1000))
#'
#' sim100 <- pred_rgev(gevp)
#' sim2 <- pred_rgev(gevp, n_years <- c(100, 1000))
#' @name pred_gev
NULL
## NULL

# ----------------------------- pred_dgev ---------------------------------

#' @rdname pred_gev
#' @export
pred_dgev <- function(ev_obj, x, n_years = 100, npy = NULL, log = FALSE) {
  #
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

#' @rdname pred_gev
#' @export
pred_pgev <- function(ev_obj, q, n_years = 100, npy = NULL,
                      lower_tail = TRUE) {
  #
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

#' @rdname pred_gev
#' @export
pred_qgev <- function(ev_obj, p, n_years = 100, npy = NULL,
                      lower_tail = TRUE) {
  #
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
    ob_fn <- function(q, ev_obj, p) {
      p_val <- pred_pgev(ev_obj = ev_obj, q = q)
      (logit(p_val) - logit(p)) ^ 2
    }
    for (j in 1:n_p) {
      qtemp <- stats::nlminb(q_mat[j, i], ob_fn, ev_obj = ev_obj, p = p[j])
      q_mat[j, i] <- qtemp$par
    }
  }
  return(q_mat)
}

# ----------------------------- pred_rgev ---------------------------------

#' @rdname pred_gev
#' @export
pred_rgev <- function(ev_obj, n_years = 100, npy = NULL) {
  #
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
    # Calculate the GEV cdf at x for each combination of (loc, scale, shape)
    # in the posterior sample, and take the mean.
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
