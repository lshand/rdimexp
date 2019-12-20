#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <cmath>

using namespace Rcpp;

//' computes empircal variogram given covariance matrix
//'
//' @param cov n x n matrix
//' @return empirical variogram n x n matrix
//' @export
// [[Rcpp::export]]
arma::mat emp_vario(arma::mat cov){
  int nobs = cov.n_rows;
  arma::mat empv;
  empv.zeros(nobs,nobs);
  for(int i=0;i<nobs;i++){
    for(int j=0;j<nobs;j++){
      empv(i,j) = cov(i,i) + cov(j,j) - 2*cov(i,j);
    }
  }
  return empv;
}

//' computes distance matrix, based on euclidean distance
//'
//' @param locs n x p matrix of location coordinates in R^p
//' @return  n x n distance matrix
//' @export
// [[Rcpp::export]]
arma::mat dist_euclid(arma::mat locs){
  int ndims = locs.n_cols;
  int nobs = locs.n_rows;

  arma::mat dist;
  dist.zeros(nobs,nobs);
  double tmp;

  for(int i=0;i<nobs;i++){
    for(int j=0;j<nobs;j++){
      tmp = 0;
      for(int k=0;k<ndims;k++){
        tmp += (locs(i,k)-locs(j,k))*(locs(i,k)-locs(j,k));
      }
      dist(i,j)  = sqrt(tmp);
    }
  }

  return dist;
}

//' Simple Kriging
//'
//' @param data n x 1 vector of observed values at locs_obs
//' @param locs_new k x p matrix of k new locations
//' @param locs_obs n x p matrix of observed locations
//' @param sig_obs n x n covariance matrix of observed data
//' @param params vector of covaraince parameter values, c(sigma2, phi)
//' @return empirical variogram, nx n matrix
//' @export
// [[Rcpp::export]]
arma::vec fitkrig(arma::vec data, arma::mat locs_new, arma::mat locs_obs,arma::mat sig_obs, arma::vec params){
  double sigma = params(0);
  double phi = params(1);
  int nobs = data.n_elem;
  int npred = locs_new.n_rows;
  arma::vec kfit;
  kfit.zeros(npred);

  locs_obs.insert_rows(0,locs_new);
  arma::mat ds_new0 = dist_euclid(locs_obs);

  arma::vec ds_new;
  ds_new.zeros(nobs);
  arma::vec kfit_vec;
  kfit_vec.zeros(nobs);
  arma::vec tmp_vec;
  tmp_vec.zeros(nobs);

  arma::vec sigdat = solve(sig_obs,data);
  double covnew;

  for(int k=0;k<npred;k++){
    ds_new = ds_new0.submat(npred,k,nobs+npred-1,k);
    kfit_vec.zeros(nobs);
    for(int i=0;i<nobs;i++){
      covnew = sigma*exp(-phi*ds_new(i));
      kfit_vec(i) = covnew * sigdat(i);
    }
    kfit(k) = sum(kfit_vec);
  }

  return kfit;
}

// [[Rcpp::Export]]
//arma::mat fitkrig_large(arma::mat data, arma::mat locs_int_all, arma::mat locs_int_clip, arma::vec params, int snipsize, double prange){
//
//  arma::mat ds = dist_euclid(locs_int_clip);
//  arma::mat Sig.s <- params(0)*exp(-params(1)*ds);
//
//   arma::mat pred;
//   pred.zeros(size(locs_int_clip));
//
//   int ncols = size(locs_int_clip,1);
//   int index0, index1l;
//
//   for(int i=0;i<ncols;i++){
//     index0 = i * (snipsize * snipsize);
//     index1 = (i+1) * snipsize * snipsize -1;
//     pred.col(i) <- fitkrig(data, locs_int_all.rows(index0,index1), locs_int_clip, Sig.s, params);
//   }
//
//   return pred;
// }

