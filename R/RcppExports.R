# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' distance to nearest pixel containing material
#'
#' @param myfield , N x 3 column matrix [X,Y,data]
#' @param matl , indiate if 0 or 1 is material
#' @return N x 1 vector of distance transforms by pixel
#' @export
dist_transform <- function(myfield, matl) {
    .Call('_rdimexp_dist_transform', PACKAGE = 'rdimexp', myfield, matl)
}

#' coarsen image
#'
#' @param myfield , N x 3 column matrix [X,Y,data]
#' @param nk , coarsen over nk x nk pixels
#' @return coarsened (N/nk x N/nk) x 3 matrix
#' @export
coarsen <- function(myfield, nk) {
    .Call('_rdimexp_coarsen', PACKAGE = 'rdimexp', myfield, nk)
}

#' grids entire image into m snips of size snipsize
#'
#' @param myfield , N x 3 column matrix [X,Y,data]
#' @param snipsize , size of snipped image (n)
#' @return n^2 x m matrix of m replicates
#' @export
make_datamat <- function(myfield, snipsize, nsmooth) {
    .Call('_rdimexp_make_datamat', PACKAGE = 'rdimexp', myfield, snipsize, nsmooth)
}

#' computes empircal variogram given covariance matrix
#'
#' @param cov n x n matrix
#' @return empirical variogram n x n matrix
#' @export
emp_vario <- function(cov) {
    .Call('_rdimexp_emp_vario', PACKAGE = 'rdimexp', cov)
}

#' computes distance matrix, based on euclidean distance
#'
#' @param locs n x p matrix of location coordinates in R^p
#' @return  n x n distance matrix
#' @export
dist_euclid <- function(locs) {
    .Call('_rdimexp_dist_euclid', PACKAGE = 'rdimexp', locs)
}

#' Simple Kriging
#'
#' @param data n x 1 vector of observed values at locs_obs
#' @param locs_new k x p matrix of k new locations
#' @param locs_obs n x p matrix of observed locations
#' @param sig_obs n x n covariance matrix of observed data
#' @param params vector of covaraince parameter values, c(sigma2, phi)
#' @return empirical variogram, nx n matrix
#' @export
fitkrig <- function(data, locs_new, locs_obs, sig_obs, params) {
    .Call('_rdimexp_fitkrig', PACKAGE = 'rdimexp', data, locs_new, locs_obs, sig_obs, params)
}

lsq_vario_fit <- function(params, empvario, ds) {
    .Call('_rdimexp_lsq_vario_fit', PACKAGE = 'rdimexp', params, empvario, ds)
}

#' profile negative log-likelihood
#'
#' @param phi
#' @param datamat n x m matrix
#' @param ds n x n matrix
#' @param prange practical range, if unsure set to max(ds)
#' @return estimated phi
#' @export
prof_nll <- function(phi, datamat, ds, prange) {
    .Call('_rdimexp_prof_nll', PACKAGE = 'rdimexp', phi, datamat, ds, prange)
}

#' optim function to estimate phi using profile neg log lik
#'
#' @param interval vector of lower and upper bounds for phi
#' @param data n x m matrix
#' @param ds n x n matrix
#' @param prange practical range, if unsure set to max(ds)
#' @return vector with estimated phi and value of obj function
#' @export
optim_nll_rcpp <- function(interval, data, ds, prange) {
    .Call('_rdimexp_optim_nll_rcpp', PACKAGE = 'rdimexp', interval, data, ds, prange)
}

#' optim function to estimate extra dimension Z
#'
#' @param init_val n+2 vector with initial values for [exp(sigma2),exp(phi),Z]
#' @param locs_obs observed coordinates n x 2 matrix
#' @param emp_vario n x n empirical variogram matrix
#' @param lambda tuning parameter lambda1
#' @return n+2 vector with estimated exponential covariance parameters and Z
#' @export
optim_z_rcpp <- function(init_val, locs_obs, emp_vario, lambda) {
    .Call('_rdimexp_optim_z_rcpp', PACKAGE = 'rdimexp', init_val, locs_obs, emp_vario, lambda)
}

#' weighted frequency of a vector
#'
#' @param x vector of values
#' @param w vector of weights
#' @return weighted contingency table
#' @export
freqweight <- function(x, w) {
    .Call('_rdimexp_freqweight', PACKAGE = 'rdimexp', x, w)
}

#' weighted frequency of a vector
#'
#' @param rvec unique distance vector
#' @param rprof radial profile
#' @param r distance matrix
#' @return empirical spectral covariance
#' @export
cov_specd <- function(rvec, rprof, r) {
    .Call('_rdimexp_cov_specd', PACKAGE = 'rdimexp', rvec, rprof, r)
}

