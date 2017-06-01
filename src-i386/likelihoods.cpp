// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// Generalized Pareto log-likelihood

// [[Rcpp::export]]
double loggp(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
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

//' Create external pointer to a C++ function for logf
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP loglik_xptr(std::string fstr) {
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "loggp")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&loggp))) ;
  else
    return(Rcpp::XPtr<funcPtr>(R_NilValue)) ;
}

// Functions for phi_to_theta, which performs transformation from phi to theta.

// GP

// [[Rcpp::export]]
Rcpp::NumericVector gp_phi_to_theta(const Rcpp::NumericVector& phi,
                                    const Rcpp::List& user_args) {
  double xm = user_args["xm"] ;
  Rcpp::NumericVector val(2);
  val[0] = phi[0] ;
  val[1] = phi[1] - phi[0] / xm ;
  return val ;
}

// A function to create external pointers to functions to evaluate
// phi_to_theta.

//' Create external pointer to a C++ function for phi_to_theta
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
