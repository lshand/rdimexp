#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <cmath>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;

arma::vec dist_euclid_dt(arma::mat locs, double x, double y){
  //int ndims = locs.n_cols;
  int nobs = locs.n_rows;

  arma::vec dist;
  dist.zeros(nobs);
  double tmp;

  //for(int i=0;i<nobs;i++){
    for(int j=0;j<nobs;j++){
      tmp = 0;
      //for(int k=0;k<ndims;k++){
      //tmp += (locs(i,k)-locs(j,k))*(locs(i,k)-locs(j,k));
      tmp += (x-locs(j,0))*(x-locs(j,0));
      tmp += (y-locs(j,1))*(y-locs(j,1));
      //}
      dist(j)  = sqrt(tmp);
    }
  //}

  return dist;
}

//' distance to nearest pixel containing material
//'
//' @param myfield , N x 3 column matrix [X,Y,data]
//' @param matl , indiate if 0 or 1 is material
//' @return N x 1 vector of distance transforms by pixel
//' @export
// [[Rcpp::export]]
arma::vec dist_transform(arma::mat myfield, int matl){
  int nn = size(myfield,0);
  arma::vec xx = myfield.col(0);
  arma::vec yy = myfield.col(1);

  arma::vec dt_vec;
  dt_vec.zeros(nn);

  //std::vector<int> z = arma::conv_to<std::vector<int>>::from(myfield.col(2));
  arma::vec z = myfield.col(2);
  arma::vec dist_vec;
  dist_vec.zeros(nn);

  //std::vector<double> dvec_tmp;
  //dvec_tmp.reserve(nn);

  Progress p(nn, true);

  //arma::vec X;
  arma::vec z_tmp, dvec_tmp;
  dvec_tmp.zeros(nn);
  z_tmp.zeros(nn);
  arma::uvec idvec, idy;
  int j;

  arma::vec ydiff;
  ydiff.zeros(sqrt(nn)*100);
  int xstart, xend;

  arma::mat myfield_tmp;

  for(int i=0; i<nn; i++){

    if(z(i)==matl){
      dt_vec(i)=0;
    }else{

      xstart = i - sqrt(nn)*50;
      xend = i + sqrt(nn)*50;
      if(xstart<0){
        xstart = 0;
      }
      if(xend>nn){
        xend = nn;
      }
      myfield_tmp = myfield.rows(xstart,xend-1);

      //ydiff =  abs(yy(i)-myfield_tmp.col(1));
      //idy = sort_index(ydiff);
      //myfield_tmp = myfield_tmp.rows(idy.subvec(0,999));

      dist_vec = dist_euclid_dt(myfield_tmp.cols(0,1), xx(i),yy(i));
      idvec = sort_index(dist_vec);
      z_tmp = myfield_tmp.col(2);

      z_tmp = z_tmp(idvec);
      dvec_tmp = dist_vec(idvec);

      j=0;
      while (z_tmp(j)==1)
      {
        j +=1;
      }

      dt_vec(i) = dvec_tmp(j);
    }

    p.increment();  // update progress
  }

  return dt_vec;
}

