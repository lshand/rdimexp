#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double lsq_vario_fit(arma::vec params, arma::mat empvario, arma::mat ds){
  double sigma = exp(params(0));
  double phi = exp(params(1));
  arma::mat vario = sigma*(1-exp(-phi*ds));
  arma::mat matsq = (empvario-vario) * (empvario-vario);
  arma::vec sumsq1 = sum(matsq,1);
  double sumsq = sum(sumsq1);
  return sumsq;
}
