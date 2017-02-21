# ================================ GEV functions ===============================

#' The Generalized Extreme Value Distribution
#'
#' Density function, distribution function, quantile function and
#' random generation for the generalized extreme value (GEV)
#' distribution with location, scale and shape parameters.
#'
#' @param x,q Numeric vectors of quantiles.
#' @param p A numeric vector of probabilities in (0,1).
#' @param loc,scale,shape Numeric vectors, of the same length.
#'   Location, scale and shape parameters.
#'   In \code{dgev}, \code{pgev} and \code{qgev} these vectors can have
#'   length > 1, but only if \code{length(x)[dgev]}, \code{length(q)[pgev]},
#'   or \code{length(p)[qgev]} has length one.
#'   Otherwise, they must have length one.
#' @param n Numeric scalar.  The number of observations to be simulated.
#' @param log A logical scalar.  If TRUE the log-density is returned.
#' @param lower_tail A logical scalar.  If TRUE (default), probabilities
#'   are P[X <= x], otherwise, P[X > x].
#' @return \code{dgev} gives the density function, \code{pgev} gives the
#'   distribution function, \code{qgev} gives the quantile function,
#'   and \code{rgev} generates random deviates.
#' @references Jenkinson, A. F. (1955) The frequency distribution of the
#'   annual maximum (or minimum) of meteorological elements.
#'   \emph{Quart. J. R. Met. Soc.}, \strong{81}, 158-171.
#'   \url{http://dx.doi.org/10.1002/qj.49708134804}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \url{http://dx.doi.org/10.1007/978-1-4471-3675-0}
#' @examples
#' dgev(-1:4, 1, 0.5, 0.8)
#' dgev(1:6, 1, 0.5, -0.2, log = TRUE)
#' dgev(1, c(1, 2), c(1, 2), c(-0.2, 0.4))
#'
#' pgev(-1:4, 1, 0.5, 0.8)
#' pgev(1:6, 1, 0.5, -0.2)
#' pgev(1, c(1, 2), c(1, 2), c(-0.2, 0.4))
#' pgev(7, c(1, 2), c(1, 2), c(-0.2, 0.4))
#' pgev(-3, c(1, 2), c(1, 2), c(-0.2, 0.4))
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

# ----------------------------- qgev ---------------------------------

#' @rdname gev
#' @export
dgev <- function (x, loc = 0, scale = 1, shape = 0,
                  log = FALSE) {
  if (min(scale) < 0) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  len_x <- length(x)
  if (len_loc != len_scale || len_scale != len_shape) {
    stop("loc, scale and shape must have the same length.")
  }
  if (len_shape != 1 & len_x != 1) {
    stop("If length(x) > 1 then loc, scale and shape must have length one.")
  }
  x <- (x - loc)/scale
  nn <- length(x)
  xx <- 1 + shape * x
  # co is a condition to ensure that calculations are only performed in
  # instances where the density is positive.
  co <- xx > 0 | is.na(xx)
  x <- x[co]
  xx <- xx[co]
  d <- numeric(nn)
  d[!co] <- -Inf
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
  if (min(scale) < 0) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  len_q <- length(q)
  if (len_loc != len_scale || len_scale != len_shape) {
    stop("loc, scale and shape must have the same length.")
  }
  if (len_shape != 1 & len_q != 1) {
    stop("If length(q) > 1 then loc, scale and shape must have length one.")
  }
  q <- (q - loc)/scale
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
  if (min(p) < 0 || max(p) > 1) {
    stop("The elements of p must all be in (0,1)")
  }
  if (min(scale) < 0) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  len_p <- length(p)
  if (len_loc != len_scale || len_scale != len_shape) {
    stop("loc, scale and shape must have the same length.")
  }
  if (len_shape != 1 & len_p != 1) {
    stop("If length(p) > 1 then loc, scale and shape must have length one.")
  }
  if (!lower_tail) {
    p <- 1 - p
  }
  xp <- -log(p)
  # If shape is close to zero then base calculation on an approximation that
  # is linear in shape.
  if (len_p == 1) {
    mult <- ifelse(abs(shape) < 1e-6, -log(xp) + shape * log(xp) ^ 2 / 2,
                   (xp ^ (-shape) - 1) / shape)
  } else {
    if(abs(shape) < 1e-6) {
      mult <- -log(xp) + shape * log(xp) ^ 2 / 2
    } else {
      mult <- (xp ^ (-shape) - 1) / shape
    }
  }
  return(loc + scale * mult)
}

# ----------------------------- qgev ---------------------------------

#' @rdname gev
#' @export
rgev <- function (n, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) {
    stop("loc, scale and shape must have length one.")
  }
  return(qgev(runif(n), loc = loc, scale = scale, shape = shape))
}

