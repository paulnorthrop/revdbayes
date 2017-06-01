// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// Generalized Pareto log-likelihood

// [[Rcpp::export]]
double cpp_gp_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  Rcpp::NumericVector gpd_data = ss["gpd_data"] ;
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

// A function to create external pointers to the functions to evaluate logf.
// See http://gallery.rcpp.org/articles/passing-cpp-function-pointers/
// If you write a new function above called new_name then add something
// like the following.
//
// else if (fstr == "new_name")
//   return(Rcpp::XPtr<funcPtr>(new funcPtr(&new_name))) ;

//' Create external pointer to a C++ log-likelihood function
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP loglik_xptr(std::string fstr) {
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gp")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&cpp_gp_loglik))) ;
  else
    return(Rcpp::XPtr<funcPtr>(R_NilValue)) ;
}

// Generalized Pareto log-prior

// [[Rcpp::export]]
double cpp_gp_mdi(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  if (x[1] < -1)
    return R_NegInf ;
  double a = pars["a"] ;
  return -log(x[0]) - a * x[1] ;
}

// [[Rcpp::export]]
double cpp_gp_flat(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  if (x[0] <= 0)
    return R_NegInf ;
  double min_xi = pars["min_xi"] ;
  double max_xi = pars["max_xi"] ;
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
SEXP prior_xptr(std::string fstr) {
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gp_mdi")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&cpp_gp_mdi))) ;
  else if (fstr == "gp_flat")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&cpp_gp_flat))) ;
  else
    return(Rcpp::XPtr<funcPtr>(R_NilValue)) ;
}

// Generalized Pareto log-posterior

// [[Rcpp::export]]
double cpp_gp_posterior(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  double loglik = cpp_gp_loglik(x, ss) ;
  double logprior = cpp_gp_mdi(x, ss) ;
  double logpost = loglik + logprior ;
  return logpost ;
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
