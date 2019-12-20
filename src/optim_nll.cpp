#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <cmath>
using namespace Rcpp;

arma::vec sigmavec(arma::mat datamat, arma::mat ds, arma::mat corsi, double prange){
  double m = size(ds,0);
  double n = size(datamat,1);
  arma::vec sigma_vec;
  sigma_vec.zeros(n);

  //for(int nobs=0;nobs<n;nobs++){
  for(int i=0; i<m; i++){
    for(int j=0; j<m; j++){
      if(ds(i,j)<=prange){
        //  sigma_vec(nobs) += 0;
        //}else{
        for(int nobs=0;nobs<n;nobs++){
          sigma_vec(nobs) += datamat(i,nobs)*datamat(j,nobs)*corsi(i,j);
        }
      }
    }
  }
  arma::vec sigma = sigma_vec/m;

  return sigma;
}

//' profile negative log-likelihood
//'
//' @param phi
//' @param datamat n x m matrix
//' @param ds n x n matrix
//' @param prange practical range, if unsure set to max(ds)
//' @return estimated phi
//' @export
// [[Rcpp::export]]
double prof_nll(double phi, arma::mat datamat, arma::mat ds, double prange){
  double m = size(ds,0);
  double n = size(datamat,1);
  arma::mat cors = exp(-phi*ds);
  arma::mat corsi = inv(cors);

  arma::vec sigma = sigmavec(datamat, ds, corsi, prange);
  double ldet = 2 * sum(log(diagvec(chol(cors))));
  double sls = 0;
  for(int i=0; i<n; i++){
    if(sigma(i)>0){
      sls +=log(sigma(i));
    }
  }
  double y = (m*sls+n*ldet+m*n)/2;
  return y;
}


//' optim function to estimate phi using profile neg log lik
//'
//' @param interval vector of lower and upper bounds for phi
//' @param data n x m matrix
//' @param ds n x n matrix
//' @param prange practical range, if unsure set to max(ds)
//' @return vector with estimated phi and value of obj function
//' @export
// [[Rcpp::export]]
arma::vec optim_nll_rcpp(arma::vec& interval, arma::mat& data, arma::mat& ds, double& prange){

  // Extract R's optim function
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optimize = stats["optimize"];

  // Call the optim function from R in C++
  Rcpp::List opt_results = optimize(Rcpp::_["f"] = Rcpp::InternalFunction(&prof_nll),
                                 //Rcpp::_["gr"]     = Rcpp::InternalFunction(&optim_gr_z_cpp),
                                 // Make sure this function is not exported!
                                 Rcpp::_["interval"]     = interval,
                                 Rcpp::_["data"] = data,
                                 Rcpp::_["ds"] = ds,
                                 Rcpp::_["prange"] = prange);

  // Extract out the estimated parameter values
  //double phi = Rcpp::as<double>(opt_results[0]);
  arma::vec out;
  out.zeros(2);
  out(0) = Rcpp::as<double>(opt_results[0]);
  out(1) = Rcpp::as<double>(opt_results[1]);

  // Return estimated value and value of obj func
  return out;
}
