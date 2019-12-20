#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <cmath>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;

//' coarsen image
//'
//' @param myfield , N x 3 column matrix [X,Y,data]
//' @param nk , coarsen over nk x nk pixels
//' @return coarsened (N/nk x N/nk) x 3 matrix
//' @export
// [[Rcpp::export]]
arma::mat coarsen(arma::mat myfield, int nk){
  arma::mat z = myfield.col(2);

  arma::vec x = myfield.col(0);
  arma::vec y = myfield.col(1);
  int ny = y.max()-y.min()+1;
  int nx = x.max()-x.min()+1;

  z.reshape(ny,nx);
  std::vector <int> xx, yy;
  int xsum = 0;
  int ysum = 0;
  int row;
  arma::mat newfield, zsub;
  newfield.zeros(ny/nk * nx/nk,3);

for(int i=0; i< nx/nk; i++){
    for(int j=0; j< ny/nk; j++){
      xx.clear();
      yy.clear();
      xsum=0;
      ysum=0;
      for(int k=0; k< nk; k++){
        xx.push_back((i+1)*nk- nk + k);
        yy.push_back((j+1)*nk- nk + k);
        xsum += xx[k];
        ysum += yy[k];
      }
      row=i*ny/nk+j;

      zsub = z.submat(yy[0],xx[0],yy[nk-1],xx[nk-1]);
      newfield(row,0) = xsum/nk;
      newfield(row,1) = ysum/nk;
      newfield(row,2) = mean(mean(zsub));
   }
  }

  return newfield;
}


//' grids entire image into m snips of size snipsize
//'
//' @param myfield , N x 3 column matrix [X,Y,data]
//' @param snipsize , size of snipped image (n)
//' @return n^2 x m matrix of m replicates
//' @export
// [[Rcpp::export]]
arma::mat make_datamat(arma::mat myfield, int snipsize, int nsmooth){
  arma::vec y = myfield.col(1);
  arma::vec x = myfield.col(0);
  int ymax = y.max()+1;
  int xmax = x.max()+1;

  int nx = xmax/snipsize;
  int ny = ymax/snipsize;
  int nn = snipsize * snipsize;

  int maxit = nx * ny;

  arma::mat datamat;
  datamat.zeros(nn/(nsmooth *nsmooth), maxit);
  Progress p(maxit, true);

  arma::mat mat_tmp_x, mat_tmp, myfield_coarse;
  mat_tmp_x.zeros(snipsize*xmax,3);
  mat_tmp.zeros(nn, 3);
  myfield_coarse.zeros(nn, 3);

  double y_tmp;
  int krow = 0;

  for(int i=0; i< nx; i++){
    mat_tmp_x = myfield.rows(i*snipsize*ymax,(i+1)*snipsize*ymax-1);

      for(int j=0; j< ny; j++){
      krow=0;
      mat_tmp.zeros(nn, 3);

      for(int k=0; k < size(mat_tmp_x,0); k++){
        y_tmp = mat_tmp_x(k,1);

        if((y_tmp < j*snipsize+snipsize) && (y_tmp > j*snipsize-1)){

          mat_tmp.row(krow) = mat_tmp_x.row(k);

          krow += 1;
        }

      }

      myfield_coarse = coarsen(mat_tmp,nsmooth);
      datamat.col(i*ny+j) = myfield_coarse.col(2);

      p.increment();  // update progress
    }

  }

  return  datamat;
}




