#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <cmath>
using namespace Rcpp;

double optim_fn_z_cpp(arma::vec params, arma::mat locs_obs, arma::mat emp_vario, const double lambda) {
  double phi = exp(params(1));
  double sig2 = exp(params(0));
  double dsz;
  double vtheta;
  double SSvariog = 0;
  double zsq = 0;
  double value = 0;

  arma::vec Z = params.subvec(2,params.n_elem-1);

  int nobs = locs_obs.n_rows;

  for(int i=0;i<nobs;i++){
    zsq += Z(i) * Z(i);
    for(int j=0;j<nobs;j++){
      dsz = 0;
      dsz += (locs_obs(i,0)-locs_obs(j,0))*(locs_obs(i,0)-locs_obs(j,0));
      dsz += (locs_obs(i,1)-locs_obs(j,1))* (locs_obs(i,1)-locs_obs(j,1));
      dsz += (Z(i)-Z(j)) * (Z(i)-Z(j));
      vtheta = sig2 * (1-exp(-phi * sqrt(dsz)));
      if(std::isnan(vtheta)){
        SSvariog += 0;
      }
      else {
        SSvariog += (emp_vario(i,j)-vtheta)*(emp_vario(i,j)-vtheta);
      }
    }

    double vnorm = sqrt(zsq);

    value = SSvariog+ lambda * vnorm;

  }

  return value;
}

arma::vec optim_gr_z_cpp(arma::vec params, arma::mat locs_obs, arma::mat emp_vario, const double lambda) {
  double phi = exp(params(1));
  double sig2 = exp(params(0));
  arma::vec Z = params.subvec(2,params.n_elem-1);

  int nobs = locs_obs.n_rows;
  int ncoef = nobs;

  //intialize derivatives;
  arma::vec dlist;
  dlist.zeros(ncoef+2);

  //calculate distances and derivatives together

  double ephi;
  arma::mat dv_z;
  dv_z.ones(nobs,nobs);
  double dz;
  double dsz;
  double tmp_dsz;
  double tmp_sig2;
  double tmp_phi;
  double tmp_z;
  double cdiff;
  double ssz = 0;

  arma::vec dz_part1;
  dz_part1.zeros(ncoef);

  for(int i=0;i<nobs;i++){
    for(int j=0;j<nobs;j++){

      tmp_dsz = 0;
      tmp_dsz += (locs_obs(i,0)-locs_obs(j,0))*(locs_obs(i,0)-locs_obs(j,0));
      tmp_dsz += (locs_obs(i,1)-locs_obs(j,1))*(locs_obs(i,1)-locs_obs(j,1));
      tmp_dsz += (Z(i)-Z(j)) * (Z(i)-Z(j));
      dsz = sqrt(tmp_dsz);
      ephi = exp(-phi * dsz);
      cdiff = (sig2*(1-ephi)-emp_vario(i,j));
      tmp_sig2 =  2*(1-ephi)*cdiff * sig2;
      tmp_phi = 2*sig2*dsz*ephi*cdiff *phi;

      if(std::isnan(ephi)){
        dlist(0) += 0; //sigma2
        dlist(1) += 0; //phi
      }
      else{
        dlist(0) += tmp_sig2;//sigma2
        dlist(1) += tmp_phi; //phi
      }

      dz = Z(i)-Z(j);
      dv_z(i,j) = 2*cdiff*sig2*phi*ephi*dz/dsz;

    } //end of j loop

    ssz += Z(i)*Z(i);

  } //end of i loop

  for(int i=0;i<nobs;i++){
    for(int j=0;j<nobs;j++){
      tmp_z = dv_z(i,j);
      if(std::isnan(tmp_z)){
        dz_part1(i) += 0;
      }
      else{
        dz_part1(i) += 2*tmp_z;
      }
    }
    dlist(i+2) = dz_part1(i) + lambda * 1/sqrt(ssz) * Z(i);
  }

  return dlist;
}

//' optim function to estimate extra dimension Z
//'
//' @param init_val n+2 vector with initial values for [exp(sigma2),exp(phi),Z]
//' @param locs_obs observed coordinates n x 2 matrix
//' @param emp_vario n x n empirical variogram matrix
//' @param lambda tuning parameter lambda1
//' @return n+2 vector with estimated exponential covariance parameters and Z
//' @export
// [[Rcpp::export]]
arma::vec optim_z_rcpp(arma::vec& init_val,
                     arma::mat& locs_obs, arma::mat& emp_vario, double& lambda){

  // Extract R's optim function
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];

  // Call the optim function from R in C++
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = init_val,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&optim_fn_z_cpp),
                                 Rcpp::_["gr"]     = Rcpp::InternalFunction(&optim_gr_z_cpp),
                                 Rcpp::_["method"] = "BFGS",
                                 Rcpp::_["control"] = "list(maxit=500)",
                                 Rcpp::_["locs_obs"] = locs_obs,
                                 Rcpp::_["emp_vario"] = emp_vario,
                                 Rcpp::_["lambda"] = lambda);

  // Extract out the estimated parameter values
  arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);

  // Return estimated values
  return out;
}
