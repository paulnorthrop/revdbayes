// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

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
  Rcpp::NumericVector zz = 1 + x[1] * sdat ;
  if (std::abs(x[1]) > 1e-6) {
    loglik = -m * log(x[0]) - (1.0 + 1.0 / x[1]) * sum(log(zz)) ;
  } else {
    double t1, t2, sdatj ;
    double total = 0;
    for(int j = 0; j < m; ++j) {
      sdatj = sdat[j] ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        total += pow(-1, i) * t1 * t2 * pow(x[1], i) / i / (i + 1) ;
      }
    }
    loglik = -m * log(x[0]) - sum_gp / x[0] - total ;
  }
  return loglik ;
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
  else
    return(Rcpp::XPtr<loglikPtr>(R_NilValue)) ;
}

// Generalized Pareto log-prior

// [[Rcpp::export]]
double cpp_gp_mdi(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  if (x[1] < -1)
    return R_NegInf ;
  double a = ppars["a"] ;
  return -log(x[0]) - a * x[1] ;
}

// [[Rcpp::export]]
double cpp_gp_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  if (x[0] <= 0)
    return R_NegInf ;
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] < min_xi)
    return R_NegInf ;
  if (x[1] > max_xi)
    return R_NegInf ;
  return -log(x[0]) ;
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
  else if (fstr == "gp_flat")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&cpp_gp_flat))) ;
  else
    return(Rcpp::XPtr<priorPtr>(R_NilValue)) ;
}

// Generalized Pareto log-posterior

// [[Rcpp::export]]
double cpp_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  // Unwrap pointer to the log-likelihood function.
  SEXP lik_ptr = pars["lik_ptr"] ;
  SEXP prior_ptr = pars["prior_ptr"] ;
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
