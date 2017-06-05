// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// Generalized Pareto log-priors

// [[Rcpp::export]]
double user_gp_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi)
    return R_NegInf ;
  return -log(x[0]) ;
}

// Generalized Extreme Value log-priors

// [[Rcpp::export]]
double user_gev_norm(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
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
double user_gev_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  return -log(x[1]) ;
}

//' Create external pointer to a C++ function for a log-prior
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP create_prior_xptr(std::string fstr) {
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& user_args) ;
  if (fstr == "gp_flat")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&user_gp_flat))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&user_gev_norm))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&user_gev_flat))) ;
  else
    return(Rcpp::XPtr<priorPtr>(R_NilValue)) ;
}

// We could create the external pointers when this file is sourced using
// this embedded R code below and/or (re)create them using the relevant
// pointer-creation functions in an R session or R package.

/*** R
  ptr_gp_flat <- create_prior_xptr("gp_flat")
  ptr_gev_norm <- create_prior_xptr("gev_norm")
  ptr_gev_flat <- create_prior_xptr("gev_flat")
*/
