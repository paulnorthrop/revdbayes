// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
double gp_mdi(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  if (x[1] < -1)
    return R_NegInf ;
  double logprior = -log(x[0]) - x[1] - 1 ;
  return logprior ;
}
