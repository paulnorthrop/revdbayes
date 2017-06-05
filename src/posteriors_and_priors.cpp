// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// Miscellaneous functions

// [[Rcpp::export]]
bool any_nonpos(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::any(x <= 0));
}

// Generalized Pareto log-likelihood

// [[Rcpp::export]]
double cpp_gp_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  double xm = ss["xm"] ;
  if (x[0] <= 0 || x[1] <= -x[0] / xm)
    return R_NegInf ;
  double loglik ;
  Rcpp::NumericVector gpd_data = ss["data"] ;
  Rcpp::NumericVector sdat = gpd_data / x[0] ;
  Rcpp::NumericVector zz = 1 + x[1] * sdat ;
  int m = ss["m"] ;
  if (std::abs(x[1]) > 1e-6) {
    loglik = -m * log(x[0]) - (1 + 1 / x[1]) * sum(log(zz)) ;
  } else {
    double sum_gp = ss["sum_gp"] ;
    double t1, t2, sdatj ;
    double total = 0.0;
    for(int j = 0; j < m; ++j) {
      sdatj = sdat[j] ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        total += pow(-1.0, i) * t1 * t2 * pow(x[1], i) / i / (i + 1) ;
      }
    }
    loglik = -m * log(x[0]) - sum_gp / x[0] - total ;
  }
  return loglik ;
}

// Generalized Extreme Value log-likelihood

// [[Rcpp::export]]
double cpp_gev_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  if (x[1] <= 0)
    return R_NegInf ;
  Rcpp::NumericVector data = ss["data"] ;
  Rcpp::NumericVector sdat = (data - x[0]) / x[1] ;
  Rcpp::NumericVector zz = 1 + x[2] * sdat ;
  if (any_nonpos(zz))
    return R_NegInf ;
  int m = ss["m"] ;
  double val = -m * log(x[1]) ;
  if (std::abs(x[2]) > 1e-6) {
    val = val - (1 + 1 / x[2]) * sum(log(zz)) - sum(pow(zz, (-1 / x[2]))) ;
  } else {
    double sum_gev = ss["sum_gev"] ;
    double t1, t2, sdatj, temp ;
    double t0 = (sum_gev - m * x[0]) / x[1] ;
    double tot = 0.0;
    double tot2 = 0.0;
    for(int j = 0; j < m; ++j) {
      sdatj = sdat[j] ;
      temp = 0.0 ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        tot += pow(-1.0, i) * t1 * t2 * pow(x[2], i) / i / (i + 1) ;
        temp += pow(-1.0, i) * pow(sdatj, i + 1) * pow(x[2], i) / (i + 1) ;
      }
      tot2 += exp(-sdatj - temp) ;
    }
    val = val - t0 - tot - tot2 ;
  }
  return val ;
}

// Order statistics log-likelihood

// [[Rcpp::export]]
double cpp_os_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  if (x[1] <= 0)
    return R_NegInf ;
  Rcpp::NumericVector data = ss["data"] ;
  Rcpp::NumericVector sdat = (data - x[0]) / x[1] ;
  Rcpp::NumericVector zz = 1 + x[2] * sdat ;
  if (any_nonpos(zz))
    return R_NegInf ;
  Rcpp::NumericVector min_data = ss["min_data"] ;
  Rcpp::NumericVector smindat = (min_data - x[0]) / x[1] ;
  Rcpp::NumericVector zz_min = 1 + x[2] * smindat ;
  int nos = ss["nos"] ;
  double val = -nos * log(x[1]) ;
  if (std::abs(x[2]) > 1e-6) {
    val = val - (1 + 1 / x[2]) * sum(log(zz)) - sum(pow(zz_min, (-1 / x[2]))) ;
  } else {
    double sum_os = ss["sum_os"] ;
    double t1, t2, sdatj, smindatj, temp ;
    double t0 = (sum_os - nos * x[0]) / x[1] ;
    double tot = 0.0;
    double tot2 = 0.0;
    for(int j = 0; j < nos; ++j) {
      sdatj = sdat[j] ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        tot += pow(-1.0, i) * t1 * t2 * pow(x[2], i) / i / (i + 1) ;
      }
    }
    int nmax = ss["nmax"] ;
    for(int j = 0; j < nmax; ++j) {
      temp = 0.0 ;
      for(int i = 1; i < 5; ++i) {
        smindatj = smindat[j] ;
        temp += pow(-1.0, i) * pow(smindatj, i + 1) * pow(x[2], i) / (i + 1) ;
      }
      tot2 += exp(-smindatj - temp) ;
    }
      val = val - t0 - tot - tot2 ;
  }
  return val ;
}

// Point process log-likelihood

// [[Rcpp::export]]
double cpp_pp_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  if (x[1] <= 0)
    return R_NegInf ;
  double thresh = ss["thresh"] ;
  double udat = (thresh - x[0]) / x[1] ;
  double zz_u = 1 + x[2] * udat ;
  if (zz_u <= 0)
    return R_NegInf ;
  Rcpp::NumericVector data = ss["data"] ;
  Rcpp::NumericVector sdat = (data - x[0]) / x[1] ;
  Rcpp::NumericVector zz = 1 + x[2] * sdat ;
  if (any_nonpos(zz))
    return R_NegInf ;
  double n_exc = ss["n_exc"] ;
  double noy = ss["noy"] ;
  double val = -n_exc * log(x[1]) ;
  if (std::abs(x[2]) > 1e-6) {
    val = val - (1 + 1 / x[2]) * sum(log(zz)) - noy * pow(zz_u, -1 / x[2]) ;
  } else {
    double sum_pp = ss["sum_pp"] ;
    double t1, t2, sdatj ;
    double t0 = (sum_pp - n_exc * x[0]) / x[1] ;
    double tot = 0.0;
    double tot2 = 0.0;
    for(int j = 0; j < n_exc; ++j) {
      sdatj = sdat[j] ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        tot += pow(-1.0, i) * t1 * t2 * pow(x[2], i) / i / (i + 1) ;
      }
    }
    for(int i = 1; i < 5; ++i) {
      tot2 += pow(-1.0, i) * pow(udat, i + 1) * pow(x[2], i) / (i + 1) ;
    }
    val = val - t0 - tot - noy * exp(-udat - tot2) ;
  }
  return val ;
}

// Generalized Pareto log-priors

// [[Rcpp::export]]
double cpp_gp_norm(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector mean = ppars["mean"] ;
  Rcpp::NumericVector icov = ppars["icov"] ;
  double c0 = log(x[0]) - mean[0] ;
  double c1 = x[1] - mean[1] ;
  double ld = icov[0] * pow(c0, 2) + 2 * icov[1] * c0 * c1 +
    icov[2]*pow(c1, 2) ;
  return (-ld / 2 - log(x[0])) ;
}

// [[Rcpp::export]]
double cpp_gp_mdi(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  double a = ppars["a"] ;
  return -log(x[0]) - a * x[1] ;
}

// [[Rcpp::export]]
double cpp_gp_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  return -log(x[0]) ;
}

// [[Rcpp::export]]
double cpp_gp_flatflat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  return 1.0 ;
}

// [[Rcpp::export]]
double cpp_gp_jeffreys(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  return -log(x[0]) - log(1.0 + x[1]) - log(1.0 + 2.0 * x[1]) / 2.0 ;
}

// [[Rcpp::export]]
double cpp_gp_beta(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector pq = ppars["pq"] ;
  return -log(x[0]) + (pq[0] - 1.0) * log(x[1] - min_xi) +
         (pq[1] - 1.0) * log(max_xi - x[1]) ;
}

// Generalized Extreme Value log-priors

// [[Rcpp::export]]
double cpp_gev_norm(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector mean = ppars["mean"] ;
  Rcpp::NumericVector icov = ppars["icov"] ;
  double c0 = x[0] - mean[0] ;
  double c1 = log(x[1]) - mean[1] ;
  double c2 = x[2] - mean[2] ;
  double ld = icov[0] * pow(c0, 2) + 2 * icov[1] * c0 * c1 +
    2 * icov[2] * c0 * c2 + icov[3] * pow(c1, 2) + 2 * icov[4] * c1 * c2 +
    icov[5] * pow(c2, 2) ;
  return (-ld / 2 - log(x[1])) ;
}

// [[Rcpp::export]]
double cpp_gev_loglognorm(const Rcpp::NumericVector& x,
                          const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector mean = ppars["mean"] ;
  Rcpp::NumericVector icov = ppars["icov"] ;
  double c0 = log(x[0]) - mean[0] ;
  double c1 = log(x[1]) - mean[1] ;
  double c2 = x[2] - mean[2] ;
  double ld = icov[0] * pow(c0, 2) + 2 * icov[1] * c0 * c1 +
    2 * icov[2] * c0 * c2 + icov[3] * pow(c1, 2) + 2 * icov[4] * c1 * c2 +
    icov[5] * pow(c2, 2) ;
  return (-ld / 2 - log(x[1]) - log(x[0])) ;
}

// [[Rcpp::export]]
double cpp_gev_mdi(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  double a = ppars["a"] ;
  return -log(x[1]) - a * x[2] ;
}

// [[Rcpp::export]]
double cpp_gev_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  return -log(x[1]) ;
}

// [[Rcpp::export]]
double cpp_gev_flatflat(const Rcpp::NumericVector& x,
                        const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  return 1.0 ;
}

// [[Rcpp::export]]
double cpp_gev_beta(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector pq = ppars["pq"] ;
  return -log(x[1]) + (pq[0] - 1.0) * log(x[2] - min_xi) +
    (pq[1] - 1.0) * log(max_xi - x[2]) ;
}

// Model-specific log-posteriors for a user-defined prior

// [[Rcpp::export]]
double gp_user_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  SEXP prior_ptr = pars["prior"] ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  priorPtr priorfun = *xpfun ;
  return cpp_gp_loglik(x, pars) + priorfun(x, pars) ;
}

// [[Rcpp::export]]
double gev_user_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  SEXP prior_ptr = pars["prior_ptr"] ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  priorPtr priorfun = *xpfun ;
  return cpp_gp_loglik(x, pars) + priorfun(x, pars) ;
}

// [[Rcpp::export]]
double os_user_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  SEXP prior_ptr = pars["prior_ptr"] ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  priorPtr priorfun = *xpfun ;
  return cpp_gp_loglik(x, pars) + priorfun(x, pars) ;
}

// [[Rcpp::export]]
double pp_user_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  SEXP prior_ptr = pars["prior_ptr"] ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  priorPtr priorfun = *xpfun ;
  return cpp_gp_loglik(x, pars) + priorfun(x, pars) ;
}

// GP Posteriors for specific in-built priors

// [[Rcpp::export]]
double gp_mdi_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_mdi(x, pars) ;
}

// [[Rcpp::export]]
double gp_norm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_norm(x, pars) ;
}

// [[Rcpp::export]]
double gp_flat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_flat(x, pars) ;
}

// [[Rcpp::export]]
double gp_flatflat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_flatflat(x, pars) ;
}

// [[Rcpp::export]]
double gp_jeffreys_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_jeffreys(x, pars) ;
}

// [[Rcpp::export]]
double gp_beta_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_beta(x, pars) ;
}

// [[Rcpp::export]]
double gev_mdi_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_mdi(x, pars) ;
}

// GEV Posteriors for specific in-built priors

// [[Rcpp::export]]
double gev_norm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_norm(x, pars) ;
}

// [[Rcpp::export]]
double gev_loglognorm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_loglognorm(x, pars) ;
}

// [[Rcpp::export]]
double gev_flat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_flat(x, pars) ;
}

// [[Rcpp::export]]
double gev_flatflat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_flatflat(x, pars) ;
}

// [[Rcpp::export]]
double gev_beta_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_beta(x, pars) ;
}

// PP Posteriors for specific in-built priors

// [[Rcpp::export]]
double pp_mdi_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_mdi(x, pars) ;
}

// [[Rcpp::export]]
double pp_norm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_norm(x, pars) ;
}

// [[Rcpp::export]]
double pp_loglognorm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_loglognorm(x, pars) ;
}

// [[Rcpp::export]]
double pp_flat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_flat(x, pars) ;
}

// [[Rcpp::export]]
double pp_flatflat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_flatflat(x, pars) ;
}

// [[Rcpp::export]]
double pp_beta_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_beta(x, pars) ;
}

// OS Posteriors for specific in-built priors

// [[Rcpp::export]]
double os_mdi_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_mdi(x, pars) ;
}

// [[Rcpp::export]]
double os_norm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_norm(x, pars) ;
}

// [[Rcpp::export]]
double os_loglognorm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_loglognorm(x, pars) ;
}

// [[Rcpp::export]]
double os_flat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_flat(x, pars) ;
}

// [[Rcpp::export]]
double os_flatflat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_flatflat(x, pars) ;
}

// [[Rcpp::export]]
double os_beta_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_beta(x, pars) ;
}

// General log-posterior, after transformation from theta to phi.

// [[Rcpp::export]]
double cpp_logpost_phi(const Rcpp::NumericVector& phi,
                       const Rcpp::List& pars, const SEXP& phi_to_theta_ptr) {
  SEXP lik_ptr = pars["lik_ptr"] ;
  SEXP prior_ptr = pars["prior_ptr"] ;
  //  SEXP phi_to_theta_ptr = phi_to_theta ;
  // Unwrap pointer to the log-likelihood function.
  typedef double (*loglikPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ss) ;
  Rcpp::XPtr<loglikPtr> xlfun(lik_ptr) ;
  loglikPtr loglikfun = *xlfun ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  priorPtr priorfun = *xpfun ;
  // Unwrap pointer to the phi_to_theta function.
  typedef Rcpp::NumericVector (*p2tPtr)(const Rcpp::NumericVector& phi,
                               const Rcpp::List& user_args) ;
  Rcpp::XPtr<p2tPtr> p2tfun(phi_to_theta_ptr) ;
  p2tPtr phi_to_theta_fun = *p2tfun ;
  Rcpp::NumericVector theta = phi_to_theta_fun(phi, pars) ;
  return loglikfun(theta, pars) + priorfun(theta, pars) ;
}

//' Create external pointer to a C++ GP log-posterior function
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP gp_logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gp_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_mdi_logpost))) ;
  else if (fstr == "gp_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_norm_logpost))) ;
  else if (fstr == "gp_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_flat_logpost))) ;
  else if (fstr == "gp_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_flatflat_logpost))) ;
  else if (fstr == "gp_jeffreys")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_jeffreys_logpost))) ;
  else if (fstr == "gp_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_beta_logpost))) ;
  else if (fstr == "gp_user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_user_logpost))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

//' Create external pointer to a C++ GEV log-posterior function
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP gev_logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gev_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_mdi_logpost))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_norm_logpost))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_loglognorm_logpost))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_flat_logpost))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_flatflat_logpost))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_beta_logpost))) ;
  else if (fstr == "gev_user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_user_logpost))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

//' Create external pointer to a C++ OS log-posterior function
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP os_logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gev_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_mdi_logpost))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_norm_logpost))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_loglognorm_logpost))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_flat_logpost))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_flatflat_logpost))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_beta_logpost))) ;
  else if (fstr == "gev_user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_user_logpost))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

//' Create external pointer to a C++ PP log-posterior function
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP pp_logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gev_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_mdi_logpost))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_norm_logpost))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_loglognorm_logpost))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_flat_logpost))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_flatflat_logpost))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_beta_logpost))) ;
  else if (fstr == "gev_user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_user_logpost))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// Generalized Pareto phi_to_theta

// [[Rcpp::export]]
Rcpp::NumericVector gp_phi_to_theta(const Rcpp::NumericVector& phi,
                                    const Rcpp::List& user_args) {
  double xm = user_args["xm"] ;
  Rcpp::NumericVector val(2);
  val[0] = phi[0] ;
  val[1] = phi[1] - phi[0] / xm ;
  return val ;
}

// A function to create external pointers to C++ phi_to_theta functions to evaluate
// phi_to_theta.

//' Create external pointer to a C++ phi_to_theta function
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP phi_to_theta_xptr(std::string fstr) {
  typedef Rcpp::NumericVector (*p2tPtr)(const Rcpp::NumericVector& phi,
                               const Rcpp::List& user_args) ;
  if (fstr == "gp")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&gp_phi_to_theta))) ;
  else
    return(Rcpp::XPtr<p2tPtr>(R_NilValue)) ;
}
