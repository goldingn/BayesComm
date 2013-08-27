#include "rtnorm.h"

using namespace Rcpp ;

SEXP rtnorm(SEXP N, SEXP MU, SEXP SI, SEXP LOW, SEXP UP){

  int n = as<int>(N);
  NumericVector mu(MU);
  NumericVector si(SI);
  NumericVector low(LOW);
  NumericVector up(UP);
  NumericVector plow = pnorm((low - mu) / si);
  NumericVector pup = pnorm((up - mu) / si);
  NumericVector ans = plow + runif(n, 0.0, 1.0) * (pup - plow);
  ans = pmax(ans, 0.0000000001);
  ans = pmin(ans, 0.999999999);
  return wrap(mu + si * qnorm(ans));
	
}

