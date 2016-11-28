# =========================== plot.evpost ===========================

#' Plot diagnostics for an evpost object
#'
#' \code{plot} method for class "evpost".  For \code{d = 1} a histogram of the
#'   simulated values is plotted with a the density function superimposed.
#'   The density is normalized crudely using the trapezium rule.  For
#'   \code{d = 2} a scatter plot of the simulated values is produced with
#'   density contours superimposed.  For \code{d > 2} pairwise plots of the
#'   simulated values are produced.
#'
#' @param x an object of class "evpost", a result of a call to
#'   \code{\link{rpost}}.
#' @param y Not used.
#' @param ... Additional arguments passed on to \code{hist}, \code{lines},
#'   \code{contour} or \code{points}.
#' @param n A numeric scalar.  Only relevant if \code{x$d = 1} or
#'   \code{x$d = 2}. The meaning depends on the value of x$d.
#' \itemize{
#'   \item {For d = 1 : n + 1 is the number of abscissae in the trapezium
#'      method used to normalize the density.}
#'   \item {For d = 2 : an n by n regular grid is used to contour the density.}
#' }
#' @param prob Numeric vector. Only relevant for d = 2.  The contour lines are
#'   drawn such that the respective probabilities that the variable lies
#'   within the contour are approximately prob.
#' @param ru_scale A logical scalar.  Should we plot data and density on the
#'   scale used in the ratio-of-uniforms algorithm (TRUE) or on the original
#'   scale (FALSE)?
#' @param rows A numeric scalar.  When \code{d} > 2 this sets the number of
#'   rows of plots.  If the user doesn't provide this then it is set
#'   internally.
#' @param xlabs,ylabs Numeric vectors.  When \code{d} > 2 these set the labels
#'   on the x and y axes respectively.  If the user doesn't provide these then
#'   the column names of the simulated data matrix to be plotted are used.
#' @param pu_only Only produce a plot relating to the posterior distribution
#'   for the threshold exceedance probability \eqn{p}. Only relevant when
#'   \code{model == "bingp"} was used in the call to \code{rpost}.
#' @param add_pu Before producing the plots add the threshold exceedance
#'   probability \eqn{p} to the parameters of the extreme value model. Only
#'   relevant when \code{model == "bingp"} was used in the call to
#'   \code{rpost}.
#' @details
#' Note that \code{suppressWarnings} is used to avoid potential benign warnings
#'   caused by passing unused graphical parameters to \code{hist} and
#'   \code{lines} via \code{...}.
#' @examples
#' # GP posterior
#' data(gom)
#' u <- stats::quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' gpg <- rpost(n = 1000, model = "gp", prior = fp, thresh = u, data = gom)
#' plot(gpg)
#'
#' @seealso \code{\link{summary.evpost}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @export
plot.evpost <- function(x, y, ..., n = ifelse(x$d == 1, 1001, 101),
                    prob = c(0.1, 0.25, 0.5, 0.75, 0.95, 0.99),
                    ru_scale = FALSE, rows = NULL, xlabs = NULL,
                    ylabs = NULL, pu_only = FALSE, add_pu = FALSE) {
  if (!inherits(x, "evpost")) {
    stop("use only with \"evpost\" objects")
  }
  if (n < 1) {
    stop("n must be no smaller than 1")
  }
  if (ru_scale) {
    plot_data <- x$sim_vals_rho
    plot_density <- x$logf_rho
  } else {
    plot_data <- x$sim_vals
    plot_density <- x$logf
  }
  #
  if (pu_only & is.null(x$bin_sim_vals)) {
    warning("pu_only = TRUE is not relevant and has been ignored",
            immediate. = TRUE)
    pu_only <- FALSE
  }
  if (add_pu & is.null(x$bin_sim_vals)) {
    warning("add_pu = TRUE is not relevant and has been ignored",
            immediate. = TRUE)
    add_pu <- FALSE
  }
  if (pu_only) {
    plot_data <- x$bin_sim_vals
    plot_density <- x$bin_logf
    x$logf_args <- x$bin_logf_args
    x$d <- 1
  }
  if (add_pu) {
    plot_data <- cbind(x$bin_sim_vals, plot_data)
    x$d <- x$d + 1
  }
  if (x$d == 1) {
    temp <- suppressWarnings(graphics::hist(plot_data, prob = TRUE,
                                            plot = FALSE))
    a <- temp$breaks[1]
    b <- temp$breaks[length(temp$breaks)]
    h <- (b-a)/n
    xx <- seq(a, b, by = h)
    density_fun <- function(z) {
      density_list <- c(list(z), x$logf_args)
      exp(do.call(plot_density, density_list))
    }
    yy <- sapply(xx, density_fun)
    # Remove any infinite, missing, or undefined values
    cond <- is.finite(yy)
    yy <- yy[cond]
    xx <- xx[cond]
    n <- length(yy) - 1
    #
    if (pu_only) {
      area <- 1
    } else {
      area <- h * (yy[1] / 2 + sum(yy[2:n]) + yy[n + 1] / 2)
    }
    yy <- yy / area
    max_y <- max(temp$density, yy)
    temp <- list(...)
    if (is.null(temp$xlab)) {
      suppressWarnings(graphics::hist(plot_data, prob = TRUE, main="",
                                      ylim = c(0, max_y), xlab = "", ...))
      if (!is.null(colnames(plot_data))) {
        graphics::title(xlab = parse(text = colnames(plot_data)[1]))
      }
    } else {
      suppressWarnings(graphics::hist(plot_data, prob = TRUE, main="",
                                      ylim = c(0, max_y), ...))
    }
    suppressWarnings(graphics::lines(xx, yy, ...))
  }
  if (x$d == 2) {
    r <- apply(plot_data, 2, range)
    xx <- seq(r[1,1], r[2,1], len = n)
    yy <- seq(r[1,2], r[2,2], len = n)
    xy <- cbind(xx, yy)
    zz <- matrix(NA, ncol = length(xx), nrow = length(yy))
    for (i in 1:length(xx)) {
      for (j in 1:length(yy)) {
        for_logf <- c(list(c(xx[i], yy[j])), x$logf_args)
        zz[i, j] <- exp(do.call(plot_density, for_logf))
      }
    }
    zz[zz == -Inf] <- NA
    dx <- diff(xx[1:2]); dy <- diff(yy[1:2]); sz <- sort(zz)
    c1 <- cumsum(sz) * dx * dy; c1 <- c1/max(c1)
    con.levs <- sapply(prob, function(x) stats::approx(c1, sz, xout = 1 - x)$y)
    #
    graphics::contour(xx, yy, zz, levels = con.levs, add = F, ann = F,
      labels = prob * 100, ...)
    graphics::points(plot_data, col = 8, ...)
    graphics::contour(xx, yy, zz, levels = con.levs, add = T, ann = T,
      labels = prob * 100, ...)
    temp <- list(...)
    if (is.null(temp$xlab)) {
      if (!is.null(colnames(plot_data))) {
        graphics::title(xlab = parse(text = colnames(plot_data)[1]))
      }
    }
    if (is.null(temp$ylab)) {
      if (!is.null(colnames(plot_data))) {
        graphics::title(ylab = parse(text = colnames(plot_data)[2]))
      }
    }
  }
  if (x$d > 2) {
    if (is.null(rows)) {
      rows <- x$d -2
    }
    cols <- ceiling(choose(x$d, 2) / rows)
    temp <- list(...)
    if (is.null(xlabs)) {
      xlabs <- colnames(plot_data)
    }
    if (is.null(ylabs)) {
      ylabs <- colnames(plot_data)
    }
    def.par <- graphics::par(no.readonly = TRUE)
    graphics::par(mfrow = c(rows, cols))
    pairwise_plots <- function(x) {
      for (i in 1:(ncol(x)-1)) {
        for (j in (i+1):ncol(x)) {
          graphics::plot(x[, i], x[, j], xlab = "", ylab = "", ...)
          graphics::title(xlab = parse(text = xlabs[i]), ylab =
                            parse(text = ylabs[j]))
        }
      }
    }
    pairwise_plots(plot_data)
    graphics::par(def.par)
  }
}

# =========================== summary.evpost ===========================

#' Summarizing an evpost object
#'
#' \code{summary} method for class "evpost"
#'
#' @param object an object of class "evpost", a result of a call to
#'   \code{\link{rpost}}.
#' @param add_pu Includes in the summary of the simulated values the threshold
#'   exceedance probability \eqn{p}. Only relevant when \code{model == "bingp"}
#'   was used in the call to \code{rpost}.
#' @param ... Additional arguments passed on to \code{print} or \code{summary}.
#' @return Prints
#' \itemize{
#'   \item {information about the ratio-of-uniforms bounding box, i.e.
#'     \code{object$box}}
#'   \item {an estimate of the probability of acceptance, i.e.
#'     \code{object$pa}}
#'   \item {a summary of the simulated values, via
#'     \code{summary(object$sim_vals)}}
#' }
#' @examples
#' # GP posterior
#' data(gom)
#' u <- stats::quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' gpg <- rpost(n = 1000, model = "gp", prior = fp, thresh = u, data = gom)
#' summary(gpg)
#' @seealso \code{\link{ru}} for descriptions of \code{object$sim_vals} and
#'   \code{object$box}.
#' @seealso \code{\link{plot.ru}} for a diagnostic plot.
#' @export
summary.evpost <- function(object, add_pu = FALSE, ...) {
  if (!inherits(object, "evpost")) {
    stop("use only with \"evpost\" objects")
  }
  cat("ru bounding box: ", "\n")
  print(object$box, ...)
  cat("\n")
  cat("estimated probability of acceptance: ", "\n")
  print(object$pa, ...)
  cat("\n")
  cat("sample summary", "\n")
  if (add_pu & is.null(object$bin_sim_vals)) {
    warning("add_pu = TRUE is not relevant and has been ignored",
            immediate. = FALSE)
    add_pu <- FALSE
  }
  if (!add_pu) {
    print(summary(object$sim_vals, ...), ...)
  } else {
    print(summary(cbind(object$bin_sim_vals, object$sim_vals), ...), ...)
  }
}

