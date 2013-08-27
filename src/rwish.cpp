#include "rwish.h"

using namespace Rcpp ;

SEXP rwish(SEXP s, SEXP DF){

	NumericVector df(DF);
	arma::mat S = as<arma::mat>(s);
	arma::mat CC = arma::chol(S);
	int n = S.n_cols;
	arma::mat a = arma::randn(n, n);
	for(int i = 0; i < n; i++)
	{
		a.diag()[i] = sqrt(as<double>(rchisq(1, df[i])));
	}
	a = arma::trimatu(a);
	a = a * CC;
	a = arma::trans(a) * a;
	return wrap(a);
	
}

