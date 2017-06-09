// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// GEV functions

//' GEV density function
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dgev_cpp(const Rcpp::NumericVector& x, const double& loc,
                             const double& scale, const double& shape) {
  if (scale <= 0.0) {
    stop("invalid scale: scale must be positive.") ;
  }
  Rcpp::NumericVector xs = (x - loc) / scale ;
  Rcpp::NumericVector d = 1 + shape * xs ;
  for(int i = 0; i < x.size(); ++i) {
    if (d[i] < 0.0) {
      d[i] = R_NegInf ;
    } else {
      if (std::abs(shape) > 1e-6) {
        Rcpp::Rcout << "shape = " << shape << std::endl;
        d[i] = -(1 + 1 / shape) * log(d[i]) - pow(d[i], -1 / shape) ;
      } else {
        d[i] = -xs[i] + shape * xs[i] * (xs[i] - 2) / 2 -
          exp(-xs[i] + shape * pow(xs[i], 2.0) / 2) ;
      }
    }
  }
  return d - log(scale) ;
}

//' GEV distribution function
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pgev_cpp(const Rcpp::NumericVector& q, const double& loc,
                             const double& scale, const double& shape) {
  if (scale <= 0.0) {
    stop("invalid scale: scale must be positive.") ;
  }
  Rcpp::NumericVector qs = (q - loc) / scale ;
  Rcpp::NumericVector p = 1 + shape * qs ;
  for(int i = 0; i < q.size(); ++i) {
    if ((std::abs(shape) > 1e-6) || (p[i] < 0.0)) {
      p[i] = exp(pow(std::max(p[i], 0.0), -1 / shape)) ;
    } else {
      p[i] = exp(-exp(-qs[i] + shape * pow(qs[i], 2.0) / 2)) ;
    }
  }
  return p ;
}
