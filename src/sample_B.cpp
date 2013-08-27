#include "sample_B.h"

using namespace Rcpp ;

SEXP sample_B(SEXP Z, SEXP x){

	arma::mat z = as<arma::mat>(Z);
	arma::mat X = as<arma::mat>(x);
	int k = X.n_cols;
	arma::mat icov = 0.1 * arma::eye(k, k);
	arma::mat XX = arma::trans(X) * X;
	arma::mat cov = arma::inv(icov + XX);
	arma::mat mn = cov * (arma::trans(X) * z);
	arma::mat ans = arma::trans(arma::randn(1, k) * chol(cov));
	return wrap(ans + mn);
	
}

