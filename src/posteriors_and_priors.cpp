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
  Rcpp::NumericVector gpd_data = ss["data"] ;
  int m = ss["m"] ;
  double xm = ss["xm"] ;
  double sum_gp = ss["sum_gp"] ;
  if (x[0] <= 0 || x[1] <= -x[0] / xm)
    return R_NegInf ;
  double loglik ;
  Rcpp::NumericVector sdat = gpd_data / x[0] ;
  Rcpp::NumericVector zz = 1.0 + x[1] * sdat ;
  if (std::abs(x[1]) > 1e-6) {
    loglik = -m * log(x[0]) - (1.0 + 1.0 / x[1]) * sum(log(zz)) ;
  } else {
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
  Rcpp::NumericVector data = ss["data"] ;
  int m = ss["m"] ;
  double sum_gev = ss["sum_gev"] ;
  if (x[1] <= 0)
    return R_NegInf ;
  Rcpp::NumericVector sdat = (data - x[0]) / x[1] ;
  Rcpp::NumericVector zz = 1.0 + x[2] * sdat ;
  if (any_nonpos(zz))
    return R_NegInf ;
  double val = -m * log(x[1]) ;
  if (std::abs(x[2]) > 1e-6) {
    val = val - (1.0 + 1.0 / x[2]) * sum(log(zz)) - sum(pow(zz, (-1.0/x[2]))) ;
  } else {
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

//' Create external pointer to a C++ log-likelihood function
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP loglik_xptr(std::string fstr) {
  typedef double (*loglikPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ss) ;
  if (fstr == "gp")
    return(Rcpp::XPtr<loglikPtr>(new loglikPtr(&cpp_gp_loglik))) ;
  else if (fstr == "gev")
    return(Rcpp::XPtr<loglikPtr>(new loglikPtr(&cpp_gev_loglik))) ;
  else
    return(Rcpp::XPtr<loglikPtr>(R_NilValue)) ;
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

//' Create external pointer to a C++ prior function
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP logprior_xptr(std::string fstr) {
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  if (fstr == "gp_mdi")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gp_mdi))) ;
  else if (fstr == "gp_norm")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gp_norm))) ;
  else if (fstr == "gp_flat")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gp_flat))) ;
  else if (fstr == "gp_flatflat")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gp_flatflat))) ;
  else if (fstr == "gp_jeffreys")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gp_jeffreys))) ;
  else if (fstr == "gp_beta")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gp_beta))) ;
  else if (fstr == "gev_mdi")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gev_mdi))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gev_norm))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gev_loglognorm))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gev_flat))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gev_flatflat))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gev_beta))) ;
  else if (fstr == "gev_mdi")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gev_mdi))) ;
  else
    return(Rcpp::XPtr<priorPtr>(R_NilValue)) ;
}

// General log-posterior

// [[Rcpp::export]]
double cpp_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  SEXP lik_ptr = pars["lik_ptr"] ;
  SEXP prior_ptr = pars["prior_ptr"] ;
  // Unwrap pointer to the log-likelihood function.
  typedef double (*loglikPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ss) ;
  Rcpp::XPtr<loglikPtr> xlfun(lik_ptr) ;
  loglikPtr loglikfun = *xlfun ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  loglikPtr priorfun = *xpfun ;
  double loglik = loglikfun(x, pars) ;
  double logprior = priorfun(x, pars) ;
  double logpost = loglik + logprior ;
  return logpost ;
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
  loglikPtr priorfun = *xpfun ;
  // Unwrap pointer to the phi_to_theta function.
  typedef Rcpp::NumericVector (*p2tPtr)(const Rcpp::NumericVector& phi,
                  const Rcpp::List& user_args) ;
  Rcpp::XPtr<p2tPtr> p2tfun(phi_to_theta_ptr) ;
  p2tPtr phi_to_theta_fun = *p2tfun ;
  Rcpp::NumericVector theta = phi_to_theta_fun(phi, pars) ;
  double loglik = loglikfun(theta, pars) ;
  double logprior = priorfun(theta, pars) ;
  double logpost = loglik + logprior ;
  return logpost ;
}

//' Create external pointer to a C++ log-posterior function
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&cpp_logpost))) ;
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
