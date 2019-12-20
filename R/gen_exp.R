#' Fit model and generate random image with stationary exponential correlation structure
#'
#' generates random image of size N x N using parameters fit to observed data of size N x N
#'
#' @param datamat - N x m matrix of m replicates of size n at observed locations
#' @param per - proportion of points n to uniformly sample across image, default p=0.05
#' @param locs_obs - coordinates of observed data, N x 2 matrix
#' @param locs_new - new indices at which to generate new image, N x 2 matrix
#' @param qthresh - quantile to threshold at
#' @param seed - random seed for generation, default NA
#' @return - prints GRF before and after thresholding (N x 4 matrix)
#' @export

fitgen_sim_stat <- function(datamat, per=0.05,locs_obs, locs_new, qthresh, seed=NA){

  if(is.na(seed)){
    seed=sample(0:1000,1)
  }

  ns=ceiling(dim(datamat)[1]*per)
  obs <- sample(1:dim(datamat)[1],ns)
  locs_obs_clip <- as.matrix(locs_obs[obs,1:2])
  datamat_clip<-datamat[obs,]
  ds_clip <- dist_euclid(locs_obs_clip)

  phi <- optim_nll_rcpp(interval=c(0,200),data=datamat_clip,ds=ds_clip,prange=max(ds_clip))[1]
  sig2 <- mean(diag(t(datamat_clip)%*%solve(exp(-phi*ds_clip),datamat_clip))/ns)
  img_new <- gen_exp(locs_obs=locs_obs_clip,locs_new=as.matrix(locs_new),params=c(sig2,phi),qthresh=qthresh ,seed=seed)
  return(list(params=c(sig2,phi),img=img_new))
}

#' Generate random image with stationary exponential correlation structure
#'
#' generates random image of size N x N using parameters fit to observed data of size N x N
#'
#' @param locs_obs - @locs_obs - coordinates of observed data, n x 2 matrix
#' @param locs_new - new indices at which to generate new image, N x 2 matrix
#' @param params - vector of exponential params, c(sigma2, phi)
#' @param qthresh - quantile to threshold at
#' @param seed - random seed for generation, default NA
#' @return - prints GRF before and after thresholding (N x 4 matrix)
#' @export

gen_exp<-function(locs_obs, locs_new, params, qthresh, seed=NA){

if(is.na(seed)){
seed=sample(0:1000,1)
}

set.seed(seed)

ns <- dim(locs_obs)[1]
ds=dist_euclid(as.matrix(locs_obs))
Sig.s <- params[1]*exp(-params[2]*ds)
data <- rnorm(ns) %*% chol(Sig.s)

ns2 <- dim(locs_new)[1]
ncl <- floor(ns2/10000)
parts <- split(1:ns2, factor(sort((1:ns2)%%ncl)))
krigparts <-lapply(1:ncl, function(x) fitkrig(data, locs_new[parts[[x]],], locs_obs, Sig.s, params))
zkrig <- unlist(krigparts)

thresh <- quantile(zkrig, qthresh)
zthresh <- ifelse(zkrig>thresh, 1,0)
newfield=as.data.frame(cbind(locs_new,zkrig,zthresh))
colnames(newfield) = c('x','y','GRF','GRF_thresh')

return(newfield)

}
