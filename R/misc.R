# =========================== process_data =====================================

process_data <- function(model, data, thresh, noy, use_noy, ros) {
  #
  # Removes missings, extracts sample summaries.
  #
  # Args:
  #   model      : character string specifying the extreme value model.
  #   data       : sample data, of a format appropriate for the model.
  #     "gp"     : vector of raw data (or, if thresh = 0, threshold excesses).
  #     "bingp"  : vector of raw data.
  #     "gev"    : vector of block maxima.
  #     "pp"     : vector of raw data.
  #     "os"     : matrix of order statistics.
  #   thresh     : extreme value threshold applied to data.
  #   use_noy    : should we use the user-supplied value of noy or the
  #                number of threshold excesses?
  #   ros        : when model == "os", the number of order statistics to
  #                retain from each row of data.
  #
  # Returns: a list containing
  #   lik_args   : basic sample summaries to add to lik_args in rpost().
  #
  lik_args <- list()
  # is.atomic(x) || is.list(x) is lik is.vector(x) but attributes are allowed
  if (model == "gp" | model == "bingp") {
    if (!(is.atomic(data) || is.list(data)) || !is.numeric(data)) {
      stop("''data'' must be a numeric vector")
    }
    nas <- is.na(data)
    data <- data[!nas]
    # Check that the threshold is not lower than the smallest observation
    if (thresh < min(data)) {
      stop("the threshold is lower than the smallest observation")
    }
    if (model == "bingp") {
      lik_args$n_raw <- length(data)              # number of raw observations
    }
    lik_args$data <- data[data > thresh] - thresh # sample threshold excesses
    if (length(lik_args$data) == 0) {
      stop("There are no data above the threshold.")
    }
    lik_args$m <- length(lik_args$data)           # sample size
    lik_args$xm <- max(lik_args$data)             # maximum threshold excess
    lik_args$sum_gp <- sum(lik_args$data)         # sum of threshold excesses
    lik_args$n_check <- lik_args$m
    return(lik_args)
  }
  if (model == "gev") {
    if (!(is.atomic(data) || is.list(data)) || !is.numeric(data)) {
      stop("''data'' must be a numeric vector")
    }
    nas <- is.na(data)
    data <- data[!nas]
    lik_args$data <- data                           # sample threshold excesses
    lik_args$m <- length(data)                      # sample size
    lik_args$x1 <- min(data)                        # minimum block maximum
    lik_args$xm <- max(data)                        # maximum block maximum
    lik_args$sum_gev <- sum(data)                   # sum of block maxima
    lik_args$n_check <- lik_args$m
    return(lik_args)
  }
  if (model == "pp") {
    if (!(is.atomic(data) || is.list(data)) || !is.numeric(data)) {
      stop("''data'' must be a numeric vector")
    }
    nas <- is.na(data)
    data <- data[!nas]
    # Check that the threshold is not lower than the smallest observation
    if (thresh < min(data)) {
      stop("the threshold is lower than the smallest observation")
    }
    lik_args$data <- data[data > thresh]    # threshold exceedances
    lik_args$n_exc <- length(lik_args$data) # number of threshold excesses
    lik_args$thresh <- thresh               # threshold
    lik_args$xm <- max(lik_args$data)       # maximum exceedance
    if (use_noy) {
      lik_args$noy <- noy                   # number of years (blocks)
    } else {
      lik_args$noy <- lik_args$n_exc        # number of years (blocks)
    }
    lik_args$sum_pp <- sum(lik_args$data)   # sum of threshold exceedances
    lik_args$m <- length(data)              # sample size
    lik_args$n_check <- lik_args$n_exc
    return(lik_args)
  }
  if (model == "os") {
    if (!is.data.frame(data) && !is.matrix(data) &&
        !(is.atomic(data) || is.list(data))) {
      stop("''data'' must be a matrix, dataframe or vector")
    }
    # Make data a matrix
    data <- as.matrix(data)
    if (!all(apply(data, 2, is.numeric))) {
      stop("''data'' must contain numeric values")
    }
    # Each row contains order statistics for a given block.
    # Remove any row that has only missing values.
    data_col <- dim(data)[2]
    nas <- !apply(is.na(data), 1, all)
    data <- data[nas, , drop = FALSE]
    # Sort each row of the data so that the first columns contains the
    # largest values, the second column the second largest, and so on.
    data <- t(apply(data, 1, sort, decreasing = TRUE, na.last = TRUE))
    # If the original data had a single column then we need to transpose the
    # vector returned from apply to create a one-column matrix.
    if (data_col == 1) {
      data <- t(data)
    }
    if (is.null(ros)) {
      ros <- ncol(data)
    } else if (ncol(data) < ros) {
      ros <- ncol(data)
      warning("data matrix has fewer than ros columns.", immediate. = TRUE)
    }
    # Use only the first ros columns
    data <- data[, 1:ros, drop = FALSE]
    # Extract a vector containing the largest value in each row.
    lik_args$max_data <- apply(data, 1, max, na.rm = TRUE)
    # Extract a vector containing the smallest value in each row.
    lik_args$min_data <- apply(data, 1, min, na.rm = TRUE)
    # Extract a vector containing all the non-missing values.
    lik_args$data <- data[!is.na(data)]
    # Calulate the number of non-missings values in each row
    lik_args$nos <- length(lik_args$data)           # sample size
    lik_args$x1 <- min(lik_args$data)               # minimum order stat
    lik_args$xm <- max(lik_args$data)               # maximum order stat
    lik_args$sum_os <- sum(lik_args$data)           # sum of all order stats
    lik_args$n_check <- nrow(data)
    return(lik_args)
  }
}

# ========================== create_ru_list ====================================

create_ru_list <- function(model, trans, rotate, min_xi, max_xi) {
  #
  # Creates a list of arguments to pass to the functions ru() or ru_rcpp()
  # in the rust package to perform ratio-of-uniforms sampling from a
  # posterior density.
  #
  # Args:
  #   model     : character string specifying the extreme value model.
  #   trans     : "none", no transformation.
  #               "BC", marginal Box-Cox transformation
  #   rotate    : if TRUE rotate posterior using Cholesky decomposition of
  #               Hessian of negated log-posterior.
  #   min_xi    : the smallest xi with a non-zero posterior density
  #   max_xi    : the largest xi with a non-zero posterior density
  #
  # Returns: a list containing the inputs model, trans, rotate and
  #   d         : the dimension of the density (number of model parameters)
  #   lower     : vector of lower bounds on the arguments of logf.
  #   upper     : vector of upper bounds on the arguments of logf.
  #   var_names : the names of the variables (posterior parameters)
  #
  if (model == "gp") {
    d <- 2L
    if (trans == "none") {
      lower <- c(0, min_xi)
      upper <- c(Inf, max_xi)
    } else if (trans == "BC") {
      lower <- c(0, 0)
      upper <- c(Inf, Inf)
    } else {
      lower <- rep(-Inf, 2)
      upper <- rep(Inf, 2)
    }
    var_names <- c("sigma[u]", "xi")
  }
  if (model == "gev" | model == "os" | model == "pp") {
    d <- 3L
    if (trans == "none") {
      lower <- c(-Inf, 0, min_xi)
      upper <- c(Inf, Inf, max_xi)
    } else if (trans == "BC") {
      lower <- c(-Inf, 0, 0)
      upper <- c(Inf, Inf, Inf)
    } else {
      lower <- rep(-Inf, 3)
      upper <- rep(Inf, 3)
    }
    var_names = c("mu","sigma", "xi")
  }
  return(list(d = d, lower = lower, upper = upper, var_names = var_names))
}

# =========================== set_which_lam ====================================

set_which_lam <- function(model) {
  #
  # Sets which_lam, the indices of the parameter vector that are to be
  # Box-Cox transformed.
  #
  # Args:
  #   model : character string specifying the extreme value model.
  #
  # Returns : a vector (of length 2)
  #
  if (model == "gp") {
    return(1:2)
  }
  if (model == "gev" | model == "pp" | model == "os") {
    return(2:3)
  }
}

# =========================== set_range_phi ====================================

set_range_phi <- function(model, phi_mid, se_phi, mult) {
  #
  # Sets min_phi and max_phi, the smallest and largest values of the
  # elements of that it is worth including in the Box-Cox grid for phi
  # given the constraints on the model parameters.
  #
  # Args:
  #   model   : character string specifying the extreme value model.
  #   phi_mid : the middle of the grid for phi.
  #   se_phi  : an estimate of the posterior standard deviation of phi.
  #   mult    : A numeric scalar.  The grid of values used to choose the
  #             Box-Cox transformation parameter lambda is based on the
  #             MAP estimate +/- mult x estimated posterior standard
  #             deviation.
  #
  # Returns : a list containing delta and
  #   min_phi : vector of smallest values of the elements of phi
  #   max_phi : vector of largest values of the elements of phi
  #   An NA indicates that no constraint is imposed.
  #
  if (model == "gp") {
    min_phi <- pmax(0, phi_mid - mult * se_phi)
    max_phi <- pmax(0, phi_mid + mult * se_phi)
    return(list(min_phi = min_phi, max_phi = max_phi))
  }
  if (model == "gev" | model == "pp" | model == "os") {
    # The first parameter is mu, which is unconstrained (hence the NAs below)
    # and is not Box-Cox transformed.
    min_phi <- pmax(c(NA, 0, 0), phi_mid - mult * se_phi, na.rm = TRUE)
    max_phi <- pmax(c(NA, 0, 0), phi_mid + mult * se_phi, na.rm = TRUE)
    return(list(min_phi = min_phi, max_phi = max_phi))
  }
}

# =========================== box_cox ===========================

box_cox <- function (x, lambda = 1, gm = 1, lambda_tol = 1e-6,
                     poly_order = 3) {
  #
  # Computes the Box-Cox transformation of a vector.
  #
  # Args:
  #   x          : A numeric vector. (Non-negative) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   gm         : A numeric scalar.  Optional scaling parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda.tol use
  #                a Taylor series expansion.
  #   poly_order : order of Taylor series polynomial in lambda used as
  #                an approximation if abs(lambda) < lambda.tol
  #
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / (lambda * gm ^ (lambda - 1))
  #
  if (abs(lambda) > lambda_tol) {
    retval <- (x ^ lambda - 1) / lambda / gm ^ (lambda - 1)
  } else if (lambda == 0) {
    retval <- log(x)
  } else if (is.infinite(x)) {
    retval <- ifelse(lambda < 0, -1 / lambda, Inf)
  } else if (x == 0) {
    retval <- ifelse(lambda > 0, -1 / lambda, -Inf)
  } else {
    i <- 0:poly_order
    retval <- sum(log(x) ^ (i+1) * lambda ^ i / factorial(i + 1))
    retval <- retval / gm ^ (lambda - 1)
  }
  return(retval)
}

# =========================== box_cox_vec ===========================

box_cox_vec <- function(x, lambda = 1, lambda_tol = 1e-6) {
  #
  # Computes the Box-Cox transformation of a vector.  If lambda is very close
  # to zero then a first order Taylor series approximation is used.
  #
  # Args:
  #   x          : A numeric vector. (Non-negative) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda.tol use
  #                a Taylor series expansion.
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / lambda
  #
  if (any(x < 0)) {
    stop("Invalid x: x must be non-negative")
  }
  max_len <- max(length(x), length(lambda))
  x <- rep_len(x, max_len)
  lambda <- rep_len(lambda, max_len)
  retval <- ifelse(abs(lambda) > lambda_tol, (x ^ lambda - 1) / lambda,
                   ifelse(lambda == 0, log(x),
                          ifelse(is.infinite(x),
                                 ifelse(lambda < 0, -1 / lambda, Inf),
                          ifelse(x == 0, ifelse(lambda > 0, -1 / lambda, -Inf),
                                 log(x) * (1 + lambda / 2)))))
  return(retval)
}

# ====================== box_cox_deriv ==========================

box_cox_deriv <- function (x, lambda = 1, lambda_tol = 1e-6,
                           poly_order = 3) {
  #
  # Computes the derivative with respect to lambda the Box-Cox
  # transformation.
  #
  # Args:
  #   x          : A numeric vector. (Positive) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda.tol use
  #                a Taylor series expansion.
  #   poly_order : order of Taylor series polynomial in lambda used as
  #                an approximation if abs(lambda) < lambda.tol
  #
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / (lambda * gm ^ (lambda - 1))
  #
  if (abs(lambda) > lambda_tol) {
    retval <- (lambda * x ^ lambda * log(x) - x ^ lambda + 1) / lambda ^ 2
  } else {
    i <- 0:poly_order
    retval <- sum(log(x) ^ (i + 2) * lambda ^ i / ((i + 2) * factorial(i)))
  }
  return(retval)
}

# ========================= check_sample_size ==================================

check_sample_size <- function(prior_name, n_check) {
  #
  # Checks that if one of the in-built improper priors is used then the sample
  # size is sufficiently large to produce a proper posterior distribution.
  # If it is not then we stop.
  #
  # Args:
  #   prior_name : A character scalar. Name of the prior.
  #   n_check    : A numeric scalar.  The sample size.
  #
  # Returns:
  #   Nothing.
  #
  if (prior_name == "gp_flat" | prior_name == "gp_flatflat") {
    if (n_check < 3) {
      stop(check_sample_size_message(prior_name, n_check))
    }
  }
  if (prior_name == "gev_flat" | prior_name == "gev_flatflat") {
    if (n_check < 4) {
      stop(check_sample_size_message(prior_name, n_check))
    }
  }
  if (prior_name == "gev_mdi" | prior_name == "gev_beta") {
    if (n_check < 2) {
      stop(check_sample_size_message(prior_name, n_check))
    }
  }
}

check_sample_size_message <- function(prior_name, n_check) {
  text1 <- "A sample size of"
  text2 <- "is not large enough to produce a proper posterior when prior"
  text3 <- "is used."
  paste(text1, n_check, text2, prior_name, text3)
}
