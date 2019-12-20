#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <cmath>
using namespace Rcpp;

//' weighted frequency of a vector
//'
//' @param x vector of values
//' @param w vector of weights
//' @return weighted contingency table
//' @export
// [[Rcpp::export]]
arma::vec freqweight(arma::vec x,arma::vec w){
  arma::vec x1 = sort(unique(x));
  arma::vec count;
  count.zeros(x1.n_elem);
  for(int i=0; i< x1.n_elem; i++){
    for(int j=0; j< x.n_elem; j++){
      if(x(j)==x1(i)){
      count(i)+= w(j);
      }
    }
  }
  return count;
}

//' weighted frequency of a vector
//'
//' @param rvec unique distance vector
//' @param rprof radial profile
//' @param r distance matrix
//' @return empirical spectral covariance
//' @export
// [[Rcpp::export]]
arma::mat cov_specd(arma::vec rvec,arma::vec rprof, arma::mat r){
  arma::mat covmat;
  covmat.zeros(size(r));

  arma::vec dd;
  double kmin;

  double n = size(r,0);
  for(int i=0; i< n; i++){
    for(int j=i; j< n; j++){
      dd = abs(r(i,j)-rvec);
      kmin = dd.index_min();
      covmat(i,j) = rprof(kmin);
      covmat(j,i) = rprof(kmin);
    }
  }

  return covmat;
}

