% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rposterior.R
\name{wbinpost}
\alias{wbinpost}
\title{Random sampling from a binomial posterior distribution, using weights}
\usage{
wbinpost(n, prior, ds_bin)
}
\arguments{
\item{n}{A numeric scalar. The size of posterior sample required.}

\item{prior}{A function to evaluate the prior, created by
\code{\link{set_bin_prior}}.
\code{prior$prior} must be \code{"bin_beta"}.}

\item{ds_bin}{A numeric list.  Sufficient statistics for inference
  about the binomial probability \eqn{p}.  Contains
\itemize{
  \item \code{sf} : a logical vector of success (\code{TRUE}) and failure
    (\code{FALSE}) indicators.
  \item \code{w} : a numeric vector of length \code{length(sf)} containing
    the values by which to multiply the observations when constructing the
    log-likelihood.
}}
}
\value{
An object (list) of class \code{"binpost"} with components

    \item{\code{bin_sim_vals}:}{An \code{n} by 1 numeric matrix of values
      simulated from the posterior for the binomial
      probability \eqn{p}}
    \item{\code{bin_logf}:}{A function returning the log-posterior for
      \eqn{p}.}
    \item{\code{bin_logf_args}:}{A list of arguments to \code{bin_logf}.}
}
\description{
Samples from the posterior distribution of the probability \eqn{p}
of a binomial distribution.  User-supplied weights are applied to each
observation when constructing the log-likelihood.
}
\details{
For \code{prior$prior == "bin_beta"} the posterior for \eqn{p}
  is a beta distribution so \code{\link[stats:Beta]{rbeta}} is used to
  sample from the posterior.
}
\examples{
u <- quantile(gom, probs = 0.65)
ds_bin <- list(sf = gom > u, w = rep(1, length(gom)))
bp <- set_bin_prior(prior = "jeffreys")
temp <- wbinpost(n = 1000, prior = bp, ds_bin = ds_bin)
graphics::hist(temp$bin_sim_vals, prob = TRUE)
}
\seealso{
\code{\link{set_bin_prior}} for setting a prior distribution
  for the binomial probability \eqn{p}.
}
