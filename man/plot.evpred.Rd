% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_evpred.R
\name{plot.evpred}
\alias{plot.evpred}
\title{Plot diagnostics for an evpred object}
\usage{
\method{plot}{evpred}(
  x,
  ...,
  leg_pos = NULL,
  leg_text = NULL,
  which_int = c("long", "short", "both")
)
}
\arguments{
\item{x}{An object of class "evpost", a result of a call to
\code{\link{rpost}}.}

\item{...}{Additional arguments passed on to \code{\link{matplot}}.}

\item{leg_pos}{A character scalar. Keyword for the position of legend.
See \code{\link{legend}}.}

\item{leg_text}{A character or expression vector.  Text for legend.
See \code{\link{legend}}.}

\item{which_int}{A character scalar.  If \code{x$type = "i"} which
intervals should be plotted?  \code{"long"} for equi-tailed intervals,
\code{"short"} for the shortest possible intervals, \code{"both"} for
both.}
}
\value{
Nothing is returned.
}
\description{
\code{plot} method for class "evpred".  Plots summarising the predictive
distribution of the largest value to be observed in N years are produced.
The plot produced depends on \code{x$type}.
If \code{x$type = "d", "p"} or \code{"q"} then
\code{\link[graphics]{matplot}} is used to produce a line plot of the
predictive density, distribution or quantile function, respectively,
with a line for each value of N in \code{x$n_years}.
If \code{x$type = "r"} then estimates of the predictive density
(from \code{\link[stats]{density}}) are plotted with a line for each N.
If \code{x$type = "i"} then lines representing estimated predictive
intervals are plotted, with the level of the interval indicated next to
the line.
}
\examples{
data(portpirie)
mat <- diag(c(10000, 10000, 100))
pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
gevp  <- rpost(n = 1000, model = "gev", prior = pn, data = portpirie)

# Predictive density function
d_gevp <- predict(gevp, type = "d", n_years = c(100, 1000))
plot(d_gevp)

# Predictive distribution function
p_gevp <- predict(gevp, type = "p", n_years = c(100, 1000))
plot(p_gevp)

# Predictive quantiles
q_gevp <- predict(gevp, type = "q", n_years = c(100, 1000))
plot(q_gevp)

# Predictive intervals
i_gevp <- predict(gevp, type = "i", n_years = c(100, 1000), hpd = TRUE)
plot(i_gevp, which_int = "both")

# Sample from predictive distribution
r_gevp <- predict(gevp, type = "r", n_years = c(100, 1000))
plot(r_gevp)
plot(r_gevp, xlim = c(4, 10))

}
\seealso{
\code{\link{predict.evpost}} for the S3 \code{predict} method
 for objects of class \code{evpost}.
}
