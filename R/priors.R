# ================================ set_prior ===============================

#' Construction of prior distributions for extreme value model parameters
#'
#' Constructs a prior distribution for use as the argument \code{prior} in
#' \code{\link{rpost}}.  The user can either specify their own prior function
#' and arguments for hyperparameters or choose from a list of in-built
#' model-specific priors.  Note that the arguments \code{model = "gev"},
#' \code{model = "pp"} and \code{model =="os"} are equivalent because
#' a prior is specified is the GEV parameterisation in each of these cases.
#' Note also that for \code{model = "pp"} the prior GEV parameterisation
#' relates to the value of \code{noy} subsequently supplied to
#' \code{\link{rpost}}.  The argument \code{model} is used for consistency
#' with \code{rpost}.
#'
#' @param prior Either
#' \itemize{
#'   \item {A function that returns the value of the prior, or}
#'   \item {A character string giving the name of the prior.
#'     See \strong{Details} for a list of priors available for each model.}
#' }
#' @param model A character string.  If \code{prior} is a character string
#'   then \code{model} gives the extreme value model to be used.  Using
#'   either \code{model = "gev"}, \code{model = "pp"} or
#'   \code{model = "os"} will result in the same (GEV) parameterisation.
#'   If \code{prior} is a function then the value of \code{model} is stored
#'   so that in the subsequent call to \code{rpost}, consistency of the
#'   prior and extreme value model parameterisations can be checked.
#' @param ... Further arguments to be passed to the user-supplied or
#'   in-built prior function.  For details of the latter see \strong{Details}.
#' @details Of the in-built named priors available in revdbayes only
#'   those specified using \code{prior = "norm"} or
#'   \code{prior = "loglognorm"} are proper.  Other proper
#'   priors are available from the package evdbayes, via the
#'   functions \code{\link[evdbayes]{prior.prob}},
#'   \code{\link[evdbayes]{prior.quant}},
#'   \code{\link[evdbayes]{prior.loglognorm}} and
#'   \code{\link[evdbayes]{prior.norm}}.  The latter is equivalent to
#'   using \code{model = "gev"} and \code{prior = "norm"} in
#'   \code{set_prior}.
#'
#'   The other in-built priors are improper, that is, the integral of the
#'   prior function over its support is not finite.  Such priors do not
#'   necessarily result in a proper posterior distribution. Northrop and
#'   Attalides (2016) consider the issue of posterior propriety in Bayesian
#'   extreme value analyses.  In most of improper priors below the prior for
#'   the scale parameter \eqn{\sigma} is taken to be \eqn{1/\sigma},
#'   i.e. a flat prior for \eqn{log \sigma}.  Here we denote the scale
#'   parameter of the GP distribution by \eqn{\sigma}, whereas we use
#'   \eqn{\sigma_u} in the revdbayes vignette.
#'
#'   For all in-bulit priors the arguments \code{min_xi} and \code{max_xi} may
#'   be supplied by the user.  The prior density is set to zero for any value
#'   of the shape parameter \eqn{\xi} that is outside
#'   (\code{min_xi}, \code{max_xi}).  This will override the default values
#'   of \code{min_xi} and \code{max_xi} in the named priors detailed above.
#'
#'   \strong{Extreme value priors.} The names of the extreme value priors
#'   set using \code{prior} and details of hyperparameters are:
#' \itemize{
#'   \item {\code{"norm"}.
#'
#'   For \code{model = "gp"}:
#'     (\eqn{log \sigma, \xi}), is bivariate normal with mean \code{mean}
#'     (a numeric vector of length 2) and covariance matrix \code{cov}
#'     (a symmetric positive definite 2 by 2 matrix).
#'
#'   For \code{model = "gev"}:
#'     (\eqn{\mu, log \sigma, \xi}), is trivariate normal with mean
#'     \code{mean} (a numeric vector of length 3) and covariance matrix
#'     \code{cov} (a symmetric positive definite 3 by 3 matrix).
#'   }
#'   \item {\code{"loglognorm"}.  For \code{model = "gev"} only:
#'     (\eqn{log \mu, log \sigma, \xi}), is trivariate normal with mean
#'     \code{mean} (a numeric vector of length 3) and covariance matrix
#'     \code{cov} (a symmetric positive definite 3 by 3 matrix).
#'   }
#'   \item {\code{"mdi"}.
#'
#'   For \code{model = "gp"}: (an extended version
#'     of) the maximal data information (MDI) prior, that is,
#'     \deqn{\pi(\sigma, \xi) = (1/ \sigma) exp[- a (\xi + 1)], for
#'     \sigma > 0, \xi >= -1, a > 0.}
#'     The value of \eqn{a} is set using the argument \code{a}.  The default
#'     value is \eqn{a = 1}, which gives the MDI prior.
#'
#'     For \code{model = "gev"}: (an extended version
#'     of) the maximal data information (MDI) prior, that is,
#'     \deqn{\pi(\mu, \sigma, \xi) = (1/ \sigma) exp[- a (\xi + 1)], for
#'     \sigma > 0, \xi >= -1, a > 0.}
#'     The value of \eqn{a} is set using the argument \code{a}.  The default
#'     value is \eqn{a = \gamma}, where \eqn{\gamma = 0.57721} is Euler's
#'     constant, which gives the MDI prior.
#'
#'     For each of these cases \eqn{\xi} must be is bounded below
#'     \emph{a priori} for the posterior to be proper
#'     (Northrop and Attalides, 2016).  An argument for the
#'     bound \eqn{\xi >= -1} is that for \eqn{\xi < -1} the GP (GEV)
#'     likelihood is unbounded above as \eqn{-\sigma/\xi}
#'     (\eqn{\mu - \sigma/\xi})) approaches the sample maximum.  In
#'     maximum likelihood estimation of GP parameters (Grimshaw, 1993)
#'     and GEV parameters a local maximum of the likelihood
#'     is sought on the region \eqn{\sigma > 0, \xi >= -1}.
#'   }
#'   \item{\code{"flat"}.
#'
#'     For \code{model = "gp"}: a flat prior for
#'     \eqn{\xi} (and for \eqn{log \sigma}):
#'     \deqn{\pi(\sigma, \xi) = (1/ \sigma), for \sigma > 0.}
#'
#'     For \code{model = "gev"}: a flat prior for
#'     \eqn{\xi} (and for \eqn{\mu} and \eqn{log \sigma}):
#'     \deqn{\pi(\mu, \sigma, \xi) = (1/ \sigma), for \sigma > 0.}
#'   }
#'   \item{\code{"flatflat"}.
#'
#'     For \code{model = "gp"}: flat priors for
#'     \eqn{\sigma} and \eqn{\xi}:
#'     \deqn{\pi(\sigma, \xi) = 1, for \sigma > 0.}
#'
#'     For \code{model = "gev"}: flat priors for \eqn{\mu}, \eqn{\sigma}
#'     and \eqn{\xi}:
#'     \deqn{\pi(\mu, \sigma, \xi) = 1, for \sigma > 0.}
#'
#'     Therefore, the posterior is proportional to the likelihood.
#'   }
#'   \item{\code{"jeffreys"}.  For \code{model = "gp"} only: the Jeffreys
#'     prior (Castellanos and Cabras, 2007):
#'     \deqn{\pi(\sigma, \xi) = 1/ [\sigma (1+\xi) \sqrt(1+2\xi)],
#'       for \sigma > 0, \xi > -1 / 2.}
#'
#'     In the GEV case the Jeffreys prior doesn't yield a proper posterior
#'     for any sample size.  See Northrop and Attalides (2016) for details.
#'   }
#'   \item{\code{"beta"}.
#'     For \code{model = "gp"}: a beta-type(p, q)
#'     prior is used for xi on the interval (\code{min_xi}, \code{max_xi}):
#'     \deqn{\pi(\sigma, \xi) = (1/\sigma) (\xi - min_xi) ^ (p-1)
#'           (max_xi - \xi) ^ (q-1), for min_xi < \xi < max_xi.}
#'
#'     For \code{model = "gev"}: similarly ...
#'     \deqn{\pi(\mu, \sigma, \xi) = (1/\sigma) (\xi - min_xi) ^ (p-1)
#'           (max_xi - \xi) ^ (q-1), for min_xi < \xi < max_xi.}
#'
#'     The argument \code{pq} is a vector containing \code{c(p,q)}.
#'     The default settings for this prior are \code{p=6,q=9} and
#'     \code{min_xi = -1/2, max_xi = 1/2}, which corresponds to the
#'     prior for \eqn{\xi} proposed in Martins and Stedinger (2000, 2001).
#'   }
#' }
#' @return A list with class \code{"evprior"}.  The first component is the
#'   input prior, i.e. either the name of the prior or a user-supplied
#'   function.  The remaining components contain the numeric values of any
#'   hyperparameters in the prior.
#' @seealso \code{\link{rpost}} for sampling from an extreme value posterior
#'   distribution.
#' @seealso \code{\link[evdbayes]{prior.prob}},
#'   \code{\link[evdbayes]{prior.quant}}, \code{\link[evdbayes]{prior.norm}}
#'   and \code{\link[evdbayes]{prior.loglognorm}} for setting a prior
#'   distribution using the evdbayes package.
#' @seealso \code{\link[evdbayes]{posterior}} for sampling from an extreme
#'   value posterior using the evdbayes package.
#' @references Castellanos, E. M. and Cabras, S. (2007) A default Bayesian
#'   procedure for the generalized Pareto distribution.
#'   \emph{Journal of Statistical Planning and Inference} \strong{137(2)},
#'   473-483. \url{http://dx.doi.org/10.1016/j.jspi.2006.01.006}.
#' @references Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
#'   for the Generalized Pareto Distribution.  \emph{Technometrics},
#'   \strong{35(2)}, 185-191. \url{http://dx.doi.org/10.2307/1269663}.
#' @references Hosking, J. R. M. and Wallis, J. R. (1987) Parameter and
#' Quantile Estimation for the Generalized Pareto Distribution. Technometrics,
#' 29(3), 339-349. \url{http://dx.doi.org/10.2307/1269343}.
#' @references Martins, E. S. and Stedinger, J. R. (2000) Generalized maximum
#'   likelihood generalized extreme value quantlie estimators for hydrologic
#'   data, \emph{Water Resources Research}, \strong{36(3)}, 737-744.
#'   \url{http://dx.doi.org/10.1029/1999WR900330}.
#' @references Martins, E. S. and Stedinger, J. R. (2001) Generalized maximum
#'   likelihood Pareto-Poisson estimators for partial duration series,
#'   \emph{Water Resources Research}, \strong{37(10)}, 2551-2557.
#'   \url{http://dx.doi.org/10.1029/2001WR000367}.
#' @references Northrop, P.J. and Attalides, N. (2016) Posterior propriety in
#'   Bayesian extreme value analyses using reference priors
#'   \emph{Statistica Sinica}, \strong{26(2)}, 721--743
#'   \url{http://dx.doi.org/10.5705/ss.2014.034}.
#' @examples
#' # Normal prior for GEV parameters (mu, log(sigma), xi).
#' mat <- diag(c(10000, 10000, 100))
#' pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
#' pn
#'
#' # Prior for GP parameters with flat prior for xi on (-1, infinity).
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' fp
#'
#' # A user-defined prior (see the vignette for details).
#' u_prior_fn <- function(x, ab) {
#'   if (x[1] <= 0 | x[2] <= -1 | x[2] >= 1) {
#'     return(-Inf)
#'   }
#'   return(-log(x[1]) + (ab[1] - 1) * log(1 + x[2]) +
#'          (ab[2] - 1) * log(1 - x[2]))
#' }
#' up <- set_prior(prior = u_prior_fn, ab = c(2, 2))
#' @export
set_prior <- function(prior = c("norm", "loglognorm", "mdi", "flat",
                                "flatflat", "jeffreys", "beta"),
                      model = c("gev", "gp", "pp", "os"), ...) {
  if (length(model) > 1) {
    warning("model not supplied: model == \"gev\" has been used.",
            immediate. = TRUE)
  }
  model <- match.arg(model)
  # If prior is a function then just return it in the required format,
  # together with any additional arguments from ....
  if (is.function(prior)) {
    temp <- list(prior = prior, ...)
    # Add trendsd to the prior list so that the prior will work with the
    # evdbayes function posterior().
    temp$trendsd <- 0
    # For the same reason we create a new prior function that has
    # trendsd as an argument.
    new_prior <- function(pars, ..., trendsd = 0) {
      return(prior(pars, ...))
    }
    temp$prior <- new_prior
    return(structure(temp, class = "evprior", model = model))
  }
  # Otherwise, call the appropriate function to set the prior with name prior.
  prior <- match.arg(prior)
  temp <- switch(model, gp = gp_prior(prior, ...), gev = gev_prior(prior, ...),
          pp = gev_prior(prior, ...), os = gev_prior(prior, ...))
  return(structure(temp, class = "evprior", model = model))
}

# ================================= GP priors =================================

gp_prior <- function(prior = c("norm", "mdi", "flat", "flatflat", "jeffreys",
                               "beta"), ...) {
  prior <- match.arg(prior)
  temp <- list(prior = paste("gp_", prior, sep=""), ...)
  # Check for unused hyperparameter names and drop them
  hpar_vec <- switch(prior, norm = c("mean", "cov"), mdi = "a_mdi",
                     flat = NULL, jeffreys = NULL, beta = "ab")
  hpar_vec <- c(hpar_vec, "min_xi", "max_xi")
  temp <- hpar_drop(temp, hpar_vec)
  # Check for problems with min_xi and/or max_xi
  if (!is.null(temp$min_xi) & !is.null(temp$max_xi)) {
    if (temp$min_xi >= temp$max_xi)
        stop("min_xi must be less than max_xi")
  }
  if (!is.null(temp$min_xi)) {
    if (prior == "mdi" & is.infinite(temp$min_xi)) {
      stop("If min_xi=-Inf then the MDI posterior is not proper: min_xi must
             be finite.")
    }
    if (prior == "jeffreys" & temp$min_xi < -1 / 2 ) {
      temp$min_xi <- -1 / 2
      warning("min_xi < -1/2 does not make sense for the Jeffreys' prior.
               min_xi = -1/2 has been used.")
    }
  }
  # Check admissibility of hyperparameters
  if (prior == "norm") {
    mean <- temp$mean
    cov <- temp$cov
    if (length(mean) != 2 | mode(mean) != "numeric")
        stop("mean must be a numeric vector of length two")
    if (!is.matrix(cov) | any(dim(cov) != 2) | mode(cov) != "numeric")
        stop("cov must be a symmetric two by two matrix")
    if (any(abs(cov - t(cov)) > .Machine$double.eps ^ 0.5))
        warning("cov may not be symmetric")
    eg <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
    if (any(eg <= 0))
        warning("cov may not be positive definite")
    icov <- solve(cov)
    icov <- icov[row(icov) >= col(icov)]
    temp$cov <- NULL
    temp <- c(temp, list(icov=icov))
  }
  if (prior == "mdi" & !is.null(temp$a_mdi)) {
    a_mdi <- temp$a_mdi
    if (length(a_mdi) != 1 | !is.numeric(a_mdi) | a_mdi <= 0)
        stop("a_mdi must be a positive numeric vector of length 1")
  }
  if (prior == "beta" & !is.null(temp$pq)) {
    pq <- temp$pq
    if (length(pq) != 2 | !is.numeric(pq) | any(pq <= 0) )
        stop("pq must be a non-negative numeric vector of length 2")
  }
  # Add trendsd to the prior list so that the prior will work with the
  # evdbayes function posterior().
  temp$trendsd <- 0
  return(temp)
}

# ------------------------------ specific GP priors -------------------------- #

gp_norm <- function(pars, mean, icov, min_xi = -Inf, max_xi = Inf,
                    trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  pars[1] <- log(pars[1])
  cpar <- pars - mean
  ld <- icov[1] * cpar[1] ^ 2 + 2 * icov[2] * cpar[1] * cpar[2]
            + icov[3] * cpar[2] ^ 2
  return(-ld / 2 - pars[1])
}

gp_mdi <- function(pars, a = 1, min_xi = -1, max_xi = Inf, trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[1]) - a * pars[2])
}

gp_flat <- function(pars, min_xi = -Inf, max_xi = Inf, trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[1]))
}

gp_flatflat <- function(pars, min_xi = -Inf, max_xi = Inf, trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  return(1)
}

gp_jeffreys <- function(pars, min_xi = -1/2, max_xi = Inf, trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[1]) - log(1 + pars[2]) - log(1 + 2 * pars[2]) / 2)
}

gp_beta <- function(pars, min_xi = -1 / 2, max_xi = 1 / 2, pq = c(6, 9),
                    trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[1]) + (pq[1] - 1) * log(pars[2] - min_xi) +
           (pq[2] - 1) * log(max_xi - pars[2]))
}

# ================================= GEV priors =================================

gev_prior <- function(prior=c("norm", "loglognorm", "mdi", "flat", "flatflat",
                              "beta"), ...) {
  prior <- match.arg(prior)
  temp <- list(prior = paste("gev_", prior, sep=""), ...)
  # Check for unused hyperparameter names and drop them
  hpar_vec <- switch(prior, norm = c("mean", "cov"), mdi = "a_mdi",
                      flat = NULL, beta = "ab")
  hpar_vec <- c(hpar_vec, "min_xi", "max_xi")
  temp <- hpar_drop(temp, hpar_vec)
  # Check for problems with min_xi and/or max_xi
  if (!is.null(temp$min_xi) & !is.null(temp$max_xi)) {
    if (temp$min_xi >= temp$max_xi)
        stop("min_xi must be less than max_xi")
  }
  if (!is.null(temp$min_xi)) {
    if (prior == "mdi" & is.infinite(temp$min_xi)) {
      stop("If min_xi=-Inf then the MDI posterior is not proper: min_xi must
            be finite.")
    }
  }
  # Check admissibility of hyperparameters
  if (prior == "norm") {
    mean <- temp$mean
    cov <- temp$cov
    if (length(mean) != 3 | mode(mean) != "numeric")
        stop("mean must be a numeric vector of length three")
    if (!is.matrix(cov) | any(dim(cov) != 3) | mode(cov) != "numeric")
        stop("cov must be a symmetric three by three matrix")
    if (any(abs(cov - t(cov)) > .Machine$double.eps ^ 0.5))
        warning("cov may not be symmetric")
    eg <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
    if (any(eg <= 0))
        warning("cov may not be positive definite")
    icov <- solve(cov)
    icov <- icov[row(icov) >= col(icov)]
    temp$cov <- NULL
    temp <- c(temp, list(icov=icov))
  }
  if (prior == "mdi" & !is.null(temp$a_mdi)) {
    a_mdi <- temp$a_mdi
    if (length(a_mdi) != 1 | !is.numeric(a_mdi) | a_mdi <= 0)
        stop("a_mdi must be a positive numeric vector of length 1")
  }
  if (prior == "beta" & !is.null(temp$pq)) {
    pq <- temp$pq
    if (length(pq) != 2 | !is.numeric(pq) | any(pq <= 0) )
        stop("pq must be a non-negative numeric vector of length 2")
  }
  # Add trendsd to the prior list so that the prior will work with the
  # evdbayes function posterior().
  temp$trendsd <- 0
  return(temp)
}

# ------------------------------ specific GEV priors ------------------------- #

gev_norm <- function(pars, mean, icov, min_xi = -Inf, max_xi = Inf,
                     trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  pars[2] <- log(pars[2])
  cpar <- pars - mean
  ld <- icov[1] * cpar[1] ^ 2 + 2 * icov[2] * cpar[1] * cpar[2]
            + 2 * icov[3] * cpar[1] * cpar[3] + icov[4] * cpar[2] ^ 2
            + 2 * icov[5] * cpar[2] * cpar[3] + icov[6] * cpar[3] ^ 2
  return(-ld / 2 - pars[2])
}

gev_loglognorm <- function(pars, mean, icov, min_xi = -Inf, max_xi = Inf,
                           trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  pars[2] <- log(pars[2])
  pars[3] <- log(pars[3])
  cpar <- pars - mean
  ld <- icov[1] * cpar[1] ^ 2 + 2 * icov[2] * cpar[1] * cpar[2]
            + 2 * icov[3] * cpar[1] * cpar[3] + icov[4] * cpar[2] ^ 2
            + 2 * icov[5] * cpar[2] * cpar[3] + icov[6] * cpar[3] ^ 2
  return(-ld / 2 - pars[2] - pars[3])
}

gev_mdi <- function(pars, a_mdi=0.5772156649015323, min_xi=-1, max_xi=Inf,
                    trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[1]) - a_mdi * pars[3])
}

gev_flat <- function(pars, min_xi = -Inf, max_xi = Inf, trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[2]))
}

gev_flatflat <- function(pars, min_xi = -Inf, max_xi = Inf, trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  return(1)
}

gev_beta <- function(pars, min_xi = -1 / 2, max_xi = 1 / 2, pq = c(6, 9)) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[2]) + (pq[1] - 1) * log(1 / 2 + pars[3]) +
           (pq[2] - 1) * log(1 / 2 - pars[3]))
}

# ================================ hpar_drop ===============================

hpar_drop <- function(x_list, hpar_vec) {
  # Check for unused hyperparameter names and drop them
  #
  # Args:
  #   x_list   : A list. A list to define a prior of class "evprior".
  #   hpar_vec : A character vector.  The names of prior hyperparameters.
  #
  # Returns:
  #   The input list x_list with any unused hyperparameters removed.
  #
  to_drop <- 1 + which(!(names(x_list[-1]) %in% hpar_vec))
  if (length(to_drop) == 1) {
    warning("This user-supplied argument is unused and has been dropped:",
            immediate. = TRUE)
  }
  if (length(to_drop) > 1) {
    warning("The following user-supplied arguments are unused and have been
            dropped:", immediate. = TRUE)
  }
  if (length(to_drop) > 0) {
    cat(names(x_list[to_drop]), "\n")
    x_list <- x_list[-to_drop]
  }
  return(x_list)
}

# ============================== set_bin_prior ================================

#' Construction of a prior distribution for a binomial probability \eqn{p}
#'
#' Constructs a prior distribution for use as the argument \code{bin_prior} in
#' \code{\link{rpost}} or in \code{\link{binpost}}.  The user can either
#' specify their own prior function and arguments for hyperparameters or choose
#' from a list of in-built priors.
#'
#' @param prior A character string giving the name of the prior for \eqn{p}.
#'   See \strong{Details} for a list of the priors available.
#' @param ... Further arguments to be passed to the user-supplied or
#'   in-built prior function.  For details of the latter see \strong{Details}.
#' @details
#'   \strong{Binomial priors.} The names of the binomial priors set using
#'   \code{bin_prior} are:
#' \itemize{
#'   \item {\code{"jeffreys"}: the \emph{Jeffreys} beta(1/2, 1/2) prior.}
#'   \item {\code{"laplace"}: the \emph{Bayes-Laplace} beta(1, 1) prior.}
#'   \item{\code{"haldane"}: the \emph{Haldane} beta(0, 0) prior.}
#'   \item{\code{"beta"}: a beta(\eqn{\alpha, \beta}) prior.  The argument
#'     \code{ab} is a vector containing \code{c}(\eqn{\alpha, \beta}).}
#'     The default is \code{ab = c(1, 1)}.
#'   \item{\code{"mdi"}: the MDI prior
#'     \eqn{\pi(p) = 1.6186 p^p (1-p)^(1-p),     0 < p < 1.}
#'   }
#' }
#' Apart from the MDI prior these are all beta distributions.
#' @return A list of class \code{"binprior"}.  The first component is the
#'   name of the input prior.  Apart from the MDI prior this will be "beta",
#'   in which case the other component of the list is a vector of length two
#'   giving the corresponding values of the beta parameters.
#' @seealso \code{\link{binpost}} for sampling from a binomial posterior
#'   distribution.
#' @examples
#' bp <- set_bin_prior(prior = "jeffreys")
#' @export
set_bin_prior <- function(prior = c("jeffreys", "laplace", "haldane", "beta",
                                    "mdi"), ...) {
  prior <- match.arg(prior)
  temp <- list(prior = paste("bin_", prior, sep=""), ...)
  # Check for unused hyperparameter names and drop them
  if (prior == "beta") {
    hpar_vec <- "ab"
  } else {
    hpar_vec <- NULL
  }
  temp <- hpar_drop(temp, hpar_vec)
  # Check admissibility of hyperparameters
  if (prior == "beta" & !is.null(temp$ab)) {
    ab <- temp$ab
    if (length(ab) != 2 | !is.numeric(ab) | any(ab <= 0) )
      stop("ab must be a non-negative numeric vector of length 2")
  }
  if (prior == "beta" & is.null(temp$ab)) {
    temp$ab <- c(1, 1)
  }
  # Jeffreys, Laplace and Haldane prior are all beta distributions.
  temp$ab <- switch(prior, jeffreys = c(1/2, 1/2), laplace = c(1, 1),
               haldane = c(0, 0), beta = temp$ab)
  if (prior %in% c("jeffreys", "laplace", "haldane")) {
    temp$prior <- paste("bin_", "beta", sep = "")
  }
  return(structure(temp, class = "binprior"))
}

