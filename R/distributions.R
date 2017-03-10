# ================================ GEV functions ===============================

#' The Generalised Extreme Value Distribution
#'
#' Density function, distribution function, quantile function and
#' random generation for the generalised extreme value (GEV)
#' distribution.
#'
#' @param x,q Numeric vectors of quantiles.
#' @param p A numeric vector of probabilities in (0,1).
#' @param loc,scale,shape Numeric vectors, of the same length.
#'   Location, scale and shape parameters.
#'   All elements of \code{scale} must be positive.
#'   In \code{dgev}, \code{pgev} and \code{qgev} these vectors can have
#'   length > 1, but only if \code{length(x)[dgev]}, \code{length(q)[pgev]},
#'   or \code{length(p)[qgev]} has length one.
#'   Then they must either have a common length or a subset of
#'   them must have a common length with the others having length one.
#'   Otherwise, they must have length one.
#' @param n Numeric scalar.  The number of observations to be simulated.
#' @param log A logical scalar.  If TRUE the log-density is returned.
#' @param lower_tail A logical scalar.  If TRUE (default), probabilities
#'   are P[X <= x], otherwise, P[X > x].
#' @details The distribution function of a GEV distribution with parameters
#'  \code{loc} = \eqn{\mu}, \code{scale} = \eqn{\sigma} (>0) and
#'  \code{shape} = \eqn{\xi} is
#'  \deqn{F(x) = exp { - [1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)} }
#'  for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  If \eqn{\xi = 0} the
#'  distribution function is defined as the limit as \eqn{\xi} tends to zero.
#'  The support of the distribution depends on \eqn{\xi}: it is
#'  \eqn{x <= \mu - \sigma / \xi} for \eqn{\xi < 0};
#'  \eqn{x >= \mu - \sigma / \xi} for \eqn{\xi > 0};
#'  and \eqn{x} is unbounded for \eqn{\xi = 0}.
#'  Note that if \eqn{\xi < -1} the GEV density function becomes infinite
#'  as \eqn{x} approaches \eqn{\mu -\sigma/\xi}.
#'
#'  See
#'  \url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
#'  for further information.
#' @return \code{dgev} gives the density function, \code{pgev} gives the
#'   distribution function, \code{qgev} gives the quantile function,
#'   and \code{rgev} generates random deviates.
#' @references Jenkinson, A. F. (1955) The frequency distribution of the
#'   annual maximum (or minimum) of meteorological elements.
#'   \emph{Quart. J. R. Met. Soc.}, \strong{81}, 158-171.
#'   Chapter 3: \url{http://dx.doi.org/10.1002/qj.49708134804}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \url{http://dx.doi.org/10.1007/978-1-4471-3675-0_3}
#' @examples
#' dgev(-1:4, 1, 0.5, 0.8)
#' dgev(1:6, 1, 0.5, -0.2, log = TRUE)
#' dgev(1, shape = c(-0.2, 0.4))
#'
#' pgev(-1:4, 1, 0.5, 0.8)
#' pgev(1:6, 1, 0.5, -0.2)
#' pgev(1, c(1, 2), c(1, 2), c(-0.2, 0.4))
#' pgev(-3, c(1, 2), c(1, 2), c(-0.2, 0.4))
#' pgev(7, 1, 1, c(-0.2, 0.4))
#'
#' qgev((1:9)/10, 2, 0.5, 0.8)
#' qgev(0.5, c(1,2), c(0.5, 1), c(-0.5, 0.5))
#'
#' p <- (1:9)/10
#' pgev(qgev(p, 1, 2, 0.8), 1, 2, 0.8)
#'
#' rgev(6, 1, 0.5, 0.8)
#' @name gev
NULL
## NULL

# ----------------------------- dgev ---------------------------------

#' @rdname gev
#' @export
dgev <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE) {
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_x <- length(x)
  if (max_len != 1 & len_x != 1) {
    stop("If length(x) > 1 then loc, scale and shape must have length 1.")
  }
  x <- (x - loc) / scale
  nn <- length(x)
  xx <- 1 + shape * x
  # co is a condition to ensure that calculations are only performed in
  # instances where the density is positive.
  co <- xx > 0 | is.na(xx)
  d <- numeric(nn)
  d[!co] <- -Inf
  d[xx == 0 & shape < -1] <- Inf
  d[xx == 0 & shape == -1] <- log(1 / scale)
  x <- x[co]
  xx <- xx[co]
  # If shape is close to zero then base calculation on approximations that
  # are linear in shape.
  if (len_x == 1) {
    shape <- shape[co]
    d[co] <- ifelse(abs(shape) < 1e-6,
                    -x + shape * x * (x - 2) / 2 - exp(-x + shape * x ^ 2 / 2),
                    -(1 + 1 / shape) * log(xx) - xx ^ (-1 / shape))
  } else {
    if(abs(shape) < 1e-6) {
      d[co] <- -x + shape * x * (x - 2) / 2 - exp(-x + shape * x ^ 2 / 2)
    } else {
      d[co] <- -(1 + 1 / shape) * log(xx) - xx ^ (-1 / shape)
    }
  }
  d <- d - log(scale)
  if (!log) {
    d <- exp(d)
  }
  return(d)
}

# ----------------------------- pgev ---------------------------------

#' @rdname gev
#' @export
pgev <- function (q, loc = 0, scale = 1, shape = 0, lower_tail = TRUE){
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_q <- length(q)
  if (max_len != 1 & len_q != 1) {
    stop("If length(q) > 1 then loc, scale and shape must have length 1.")
  }
  q <- (q - loc) / scale
  nn <- length(q)
  qq <- 1 + shape * q
  # co is a condition to ensure that calculations are only performed in
  # instances where the density is positive.
  co <- qq > 0 | is.na(qq)
  q <- q[co]
  qq <- qq[co]
  p <- numeric(nn)
  # If shape is close to zero then base calculation on an approximation that
  # is linear in shape.
  if (len_q == 1) {
    p[!co] <- ifelse(shape[!co] > 0, 0, 1)
    shape <- shape[co]
    p[co] <- ifelse(abs(shape) < 1e-6, exp(-exp(-q + shape * q ^ 2 / 2)),
                exp(-pmax(1 + shape * q, 0) ^ (-1 / shape)))
  } else {
    p[!co] <- ifelse(shape > 0, 0, 1)
    if(abs(shape) < 1e-6) {
      p[co] <- exp(-exp(-q + shape * q ^ 2 / 2))
    } else {
      p[co] <- exp(-pmax(1 + shape * q, 0) ^ (-1 / shape))
    }
  }
  if (!lower_tail) {
    p <- 1 - p
  }
  return(p)
}

# ----------------------------- qgev ---------------------------------

#' @rdname gev
#' @export
qgev <- function (p, loc = 0, scale = 1, shape = 0, lower_tail = TRUE) {
  if (any(p <= 0) || any(p >= 1)) {
    stop("The elements of p must all be in (0,1)")
  }
  if (min(scale) < 0) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_p <- length(p)
  if (max_len != 1 & len_p != 1) {
    stop("If length(p) > 1 then loc, scale and shape must have length 1.")
  }
  if (!lower_tail) {
    p <- 1 - p
  }
  xp <- -log(p)
  # If shape is close to zero then base calculation on an approximation that
  # is linear in shape.
  mult <- box_cox_vec(x = xp, lambda = -shape, poly_order = 1)
  return(loc - scale * mult)
}

# ----------------------------- rgev ---------------------------------

#' @rdname gev
#' @export
rgev <- function (n, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) {
    stop("loc, scale and shape must have length one.")
  }
  return(qgev(stats::runif(n), loc = loc, scale = scale, shape = shape))
}

# --------------------------- rgev_vec -------------------------------

# Vectorized version of rgev.

rgev_vec <- Vectorize(rgev, vectorize.args = c("n", "loc", "scale", "shape"))

# ================================ GP functions ===============================

#' The Generalised Pareto Distribution
#'
#' Density function, distribution function, quantile function and
#' random generation for the generalised Pareto (GP) distribution.
#'
#' @param x,q Numeric vectors of quantiles.  All elements of \code{x}
#'   and \code{q} must be non-negative.
#' @param p A numeric vector of probabilities.
#'   If \code{lower_tail = TRUE} then \code{p} must be in [0, 1).
#'   If \code{lower_tail = FALSE} then \code{p} must be in (0, 1].
#' @param loc,scale,shape Numeric vectors.
#'   Location, scale and shape parameters.
#'   All elements of \code{scale} must be positive.
#'   In \code{dgp}, \code{pgp} and \code{qgp} these vectors can have
#'   length > 1, but only if \code{length(x)[dgp]}, \code{length(q)[pgp]},
#'   or \code{length(p)[qgp]} has length one.
#'   Then they must either have a common length or a subset of
#'   them must have a common length with the others having length one.
#'   Otherwise, they must have length one.
#' @param n Numeric scalar.  The number of observations to be simulated.
#' @param log A logical scalar.  If TRUE the log-density is returned.
#' @param lower_tail A logical scalar.  If TRUE (default), probabilities
#'   are P[X <= x], otherwise, P[X > x].
#' @details The distribution function of a GP distribution with parameters
#'  \code{location} = \eqn{\mu}, \code{scale} = \eqn{\sigma} (>0) and
#'  \code{shape} = \eqn{\xi} is
#'  \deqn{F(x) = 1 - [1 + \xi (x - \mu) / \sigma] ^ (-1/\xi) }
#'  for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  If \eqn{\xi = 0} the
#'  distribution function is defined as the limit as \eqn{\xi} tends to zero.
#'  The support of the distribution depends on \eqn{\xi}: it is
#'  \eqn{x >= \mu} for \eqn{\xi >= 0};
#'  and \eqn{0 <= x <= \mu - \sigma / \xi} for \eqn{\xi < 0}.  Note that
#'  if \eqn{\xi < -1} the GP density function becomes infinite as \eqn{x}
#'  approaches \eqn{\mu - \sigma/\xi}.
#'
#'  See
#'  \url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#'   for further information.
#' @return \code{dgp} gives the density function, \code{pgp} gives the
#'   distribution function, \code{qgp} gives the quantile function,
#'   and \code{rgp} generates random deviates.
#' @references Pickands, J. (1975) Statistical inference using extreme
#'   order statistics. \emph{Annals of Statistics}, \strong{3}, 119-131.
#'   \url{http://dx.doi.org/10.1214/aos/1176343003}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   Chapter 4: \url{http://dx.doi.org/10.1007/978-1-4471-3675-0_4}
#' @examples
#' dgp(0:4, scale = 0.5, shape = 0.8)
#' dgp(1:6, scale = 0.5, shape = -0.2, log = TRUE)
#' dgp(1, scale = 1, shape = c(-0.2, 0.4))
#'
#' pgp(0:4, scale = 0.5, shape = 0.8)
#' pgp(1:6, scale = 0.5, shape = -0.2)
#' pgp(1, scale = c(1, 2), shape = c(-0.2, 0.4))
#' pgp(7, scale = 1, shape = c(-0.2, 0.4))
#'
#' qgp((0:9)/10, scale = 0.5, shape = 0.8)
#' qgp(0.5, scale = c(0.5, 1), shape = c(-0.5, 0.5))
#'
#' p <- (1:9)/10
#' pgp(qgp(p, scale = 2, shape = 0.8), scale = 2, shape = 0.8)
#'
#' rgp(6, scale = 0.5, shape = 0.8)
#' @name gp
NULL
## NULL

# ----------------------------- dgp ---------------------------------

#' @rdname gp
#' @export
dgp <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE) {
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_x <- length(x)
  if (max_len != 1 & len_x != 1) {
    stop("If length(x) > 1 then scale and shape must have length 1.")
  }
  x <- (x - loc) / scale
  nn <- length(x)
  xx <- 1 + shape * x
  # co is a condition to ensure that calculations are only performed in
  # instances where the density is positive.
  co <- (x >= 0 & xx > 0) | is.na(xx)
  d <- numeric(nn)
  d[!co] <- -Inf
  d[xx == 0 & shape < -1] <- Inf
  d[xx == 0 & shape == -1] <- log(1 / scale)
  x <- x[co]
  xx <- xx[co]
  # If shape is close to zero then base calculation on approximations that
  # are linear in shape.
  if (len_x == 1) {
    shape <- shape[co]
    d[co] <- ifelse(abs(shape) < 1e-6, -x + shape * x * (x - 2) / 2,
                    -(1 + 1 / shape) * log(xx))
  } else {
    if(abs(shape) < 1e-6) {
      d[co] <- -x + shape * x * (x - 2) / 2
    } else {
      d[co] <- -(1 + 1 / shape) * log(xx)
    }
  }
  d <- d - log(scale)
  if (!log) {
    d <- exp(d)
  }
  return(d)
}

# ----------------------------- pgp ---------------------------------

#' @rdname gp
#' @export
pgp <- function (q, loc = 0, scale = 1, shape = 0, lower_tail = TRUE){
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_q <- length(q)
  if (max_len != 1 & len_q != 1) {
    stop("If length(q) > 1 then scale and shape must have length 1.")
  }
  q <- pmax(q - loc, 0) / scale
  nn <- length(q)
  qq <- 1 + shape * q
  # co is a condition to ensure that calculations are only performed in
  # instances where the density is positive.
  co <- qq > 0 | is.na(qq)
  q <- q[co]
  qq <- qq[co]
  p <- numeric(nn)
  p[!co] <- 1
  # If shape is close to zero then base calculation on an approximation that
  # is linear in shape.
  if (len_q == 1) {
    shape <- shape[co]
    p[co] <- ifelse(abs(shape) < 1e-6, 1 - exp(-q + shape * q ^ 2 / 2),
                    1 - pmax(1 + shape * q, 0) ^ (-1 / shape))
  } else {
    if(abs(shape) < 1e-6) {
      p[co] <- 1 - exp(-q + shape * q ^ 2 / 2)
    } else {
      p[co] <- 1 - pmax(1 + shape * q, 0) ^ (-1 / shape)
    }
  }
  if (!lower_tail) {
    p <- 1 - p
  }
  return(p)
}

# ----------------------------- qgp ---------------------------------

#' @rdname gp
#' @export
qgp <- function (p, loc = 0, scale = 1, shape = 0, lower_tail = TRUE) {
  if (lower_tail) {
    if (any(p < 0) || any(p >= 1)) {
      stop("The elements of p must all be in [0,1)")
    }
  }
  if (!lower_tail) {
    if (any(p <= 0) || any(p > 1)) {
      stop("The elements of p must all be in (0,1]")
    }
  }
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_p <- length(p)
  if (max_len != 1 & len_p != 1) {
    stop("If length(p) > 1 then scale and shape must have length 1.")
  }
  if (!lower_tail) {
    p <- 1 - p
  }
  # If shape is close to zero then base calculation on an approximation that
  # is linear in shape.
  mult <- box_cox_vec(x = 1 - p, lambda = -shape, poly_order = 1)
  return(loc - scale * mult)
}

# --------------------------- qgp_vec -------------------------------

# Vectorized version of qgp.

  qgp_vec <- Vectorize(qgp, vectorize.args = c("p", "loc", "scale", "shape"))

# ----------------------------- rgp ---------------------------------

#' @rdname gp
#' @export
rgp <- function (n, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) {
    stop("loc, scale and shape must have length one.")
  }
  return(qgp(stats::runif(n), loc = loc, scale = scale, shape = shape))
}

# --------------------------- rgp_vec -------------------------------

# Vectorized version of rgp.

rgp_vec <- Vectorize(rgp, vectorize.args = c("n", "scale", "shape"))

# =========================== binomial-GP functions ============================

# --------------------------- dbingp ---------------------------------

dbingp <- function(x, p_u = 0.5 , loc = 0, scale = 1, shape = 0, log = FALSE) {
  #
  # Binomial-GP `density' function.
  #`
  # Args:
  #   x     : Numeric vector of quantiles.  No element of x can be < loc.
  #   p_u   : Numeric vector of threshold exceedance probabilities in (0,1).
  #   loc   : Numeric vector of GP location parameters: usually the threshold.
  #   scale : Numeric vector of GP scale parameters.
  #   shape : Numeric vector of GP shape parameters.
  #   log   : A logical scalar.  If TRUE the log-density is returned.
  #
  #   All elements of scale must be positive. These vectors can have
  #   length > 1, but only if x has length one.  Then they must either have
  #   a common length or a subset of them must have a common length with the
  #   others having length one. Otherwise, they must have length one.
  #
  if (any(x < loc)) {
    stop("Invalid x: no element of  can be less than loc.")
  }
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  if (any(p_u <= 0) || any(p_u >= 1)) {
    stop("invalid p_u: p_u must be in (0,1).")
  }
  len_p_u <- length(p_u)
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_p_u, len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("p_u, loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  p_u <- rep(p_u, length.out = max_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_x <- length(x)
  if (max_len != 1 & len_x != 1) {
    stop("If length(x) > 1 then scale and shape must have length 1.")
  }
  d <- dgp(x = x, loc = loc, scale = scale, shape = shape, log = TRUE) +
    log(p_u)
  if (!log) {
    d <- exp(d)
  }
  return(d)
}

# --------------------------- pbingp ---------------------------------

pbingp <- function(q, p_u = 0.5 , loc = 0, scale = 1, shape = 0,
                   lower_tail = TRUE) {
  #
  # Binomial-GP distribution function.
  #`
  # Args:
  #   q          : Numeric vector of quantiles.  No element of q can be < loc.
  #   p_u        : Numeric vector of threshold exceedance probabilities in
  #                (0,1).
  #   loc        : Numeric vector of GP location parameters: usually the
  #               threshold.
  #   scale      : Numeric vector of GP scale parameters.
  #   shape      : Numeric vector of GP shape parameters.
  #   lower_tail : A logical scalar.  If TRUE (default), probabilities
  #                are P[X <= x], otherwise, P[X > x].
  #
  #   All elements of scale must be positive. These vectors can have
  #   length > 1, but only if q has length one.  Then they must either have
  #   a common length or a subset of them must have a common length with the
  #   others having length one. Otherwise, they must have length one.
  #
  if (any(q < loc)) {
    stop("Invalid q: no element of q can be less than loc.")
  }
  if (min(scale) < 0) {
    stop("invalid scale: scale must be positive.")
  }
  if (any(p_u <= 0) || any(p_u >= 1)) {
    stop("invalid p_u: p_u must be in (0,1).")
  }
  len_p_u <- length(p_u)
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_p_u, len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("p_u, loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  p_u <- rep(p_u, length.out = max_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_q <- length(q)
  if (max_len != 1 & len_q != 1) {
    stop("If length(q) > 1 then scale and shape must have length 1.")
  }
  p <- 1 - p_u * pgp(q = q, loc = loc, scale = scale, shape = shape,
                     lower_tail = FALSE)
  if (!lower_tail) {
    p <- 1 - p
  }
  return(p)
}

# --------------------------- qbingp ---------------------------------

qbingp <- function(p, p_u = 0.5, loc = 0, scale = 1, shape = 0,
                   lower_tail = TRUE) {
  #
  # Binomial-GP quantiles.
  #`
  # Args:
  #   p          : Numeric vector of probabilities in (0,1).
  #   p_u        : Numeric vector of threshold exceedance probabilities in
  #                (0,1).
  #   loc        : Numeric vector of GP location parameters: usually the
  #               threshold.
  #   scale      : Numeric vector of GP scale parameters.
  #   shape      : Numeric vector of GP shape parameters.
  #   lower_tail : A logical scalar.  If TRUE (default), probabilities
  #                are P[X <= x], otherwise, P[X > x].
  #
  #   All elements of scale must be positive. These vectors can have
  #   length > 1, but only if p has length one.  Then they must either have
  #   a common length or a subset of them must have a common length with the
  #   others having length one. Otherwise, they must have length one.
  #
  if (any(p <= 0) || any(p >= 1)) {
    stop("invalid p: p must be in (0,1).")
  }
  if (!lower_tail) {
    p <- 1 - p
  }
  if (any(p < 1 - p_u)) {
    if (lower_tail) {
      stop("Invalid p: no element of p can be less than 1 - p_u.")
    }
    if (!lower_tail) {
      stop("Invalid p: no element of p can be greater than than p_u.")
    }
  }
  if (min(scale) < 0) {
    stop("invalid scale: scale must be positive.")
  }
  if (any(p_u <= 0) || any(p_u >= 1)) {
    stop("invalid p_u: p_u must be in (0,1).")
  }
  len_p_u <- length(p_u)
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_p_u, len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("p_u, loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  p_u <- rep(p_u, length.out = max_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_p <- length(p)
  if (max_len != 1 & len_p != 1) {
    stop("If length(p) > 1 then p_u, loc, scale and shape must have length 1.")
  }
  pnew <- 1 - (1 - p) / p_u
  q <- qgp_vec(p = pnew, loc = loc, scale = scale, shape = shape)
  return(q)
}
