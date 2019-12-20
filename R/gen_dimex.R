#' Fit and Generate random image with nonstationary correlations structure
#'
#' generates random image of size N x N using parameters fit to observed data of size N x N
#'
#' @param locs_obs - coordinates of observed data, n x 2 matrix
#' @param locs_new - new indices at which to generate new image, N x 2 matrix
#' @param datamat - N x m matrix of m replicates of size n at observed locations
#' @param per - proportion of points n to uniformly sample across image, default p=0.05
#' @param lambdas - list specify tuning parameters lambda1 and lambda2
#' @param qthresh - quantile to threshold at
#' @param seed - random seed for generation, default NA
#' @return - prints GRF before and after thresholding plus estimated Z (N x 5 matrix)
#' @export

fitgen_sim <- function(locs_obs, locs_new, datamat,per=0.05,lambdas, qthresh,seed){
  fit_new <- fit_dimex(locs_obs=locs_obs,datamat=datamat,per=per,lambdas = lambdas, seed=seed, extra=FALSE, progress=FALSE)
  img_new <- gen_dimex(fit_new$locs_obs,locs_new=locs_new,lambda2=lambdas$lambda2,theta=fit_new$opt_params,qthresh=qthresh ,seed=seed, par=FALSE)
  return(list(params=fit_new$opt_params,img=img_new))
}

#' Generate random image with nonstationary correlations structure
#'
#' generates random image of size N x N using parameters fit to observed data of size N x N
#'
#' @param locs_obs - coordinates of observed data, n x 2 matrix
#' @param locs_new - new indices at which to generate new image, N x 2 matrix
#' @param lambda2 - tuning parameter lambda2.
#' @param theta - [sigma2, phi, Z] params from fit_dimex()
#' @param qthresh - quantile to threshold at
#' @param seed - random seed for generation, default NA
#' @param par - true/false if to parallelize, default=FALSE. should not be TRUE if generating multiple GRFs in parallel already. should=TRUE if generating one image at a time.
#' @return - prints GRF before and after thresholding plus estimated Z (N x 5 matrix)
#' @export

gen_dimex<-function(locs_obs, locs_new, lambda2, theta, qthresh, seed=NA, par=FALSE){

if(is.na(seed)){
seed=sample(0:1000)
}

set.seed(seed)

ns <- dim(locs_obs)[1]
params = theta[1:2]

locs_int = cbind(locs_obs,theta[-(1:2)])
fit.s = fields::Tps(x=locs_int[,-3], Y=locs_int[,3],lambda=lambda2, give.warnings=FALSE)
locs_int[,3] = predict(fit.s,locs_int[,-3], locs=0)
ds=dist_euclid(as.matrix(locs_int))

zall <- predict(fit.s,locs_new)
locs_new_all <- as.matrix(cbind(locs_new,zall))

Sig.s <- params[1]*exp(-params[2]*ds)
data <- rnorm(ns) %*% chol(Sig.s)

ns2 <- dim(locs_new_all)[1]
ncl <- floor(ns2/10000)
parts <- split(1:ns2, factor(sort((1:ns2)%%ncl)))
if(par){
krigparts <-future_lapply(1:ncl, function(x) fitkrig(data, locs_new_all[parts[[x]],], locs_int, Sig.s, params))
}else{
krigparts <-lapply(1:ncl, function(x) fitkrig(data, locs_new_all[parts[[x]],], locs_int, Sig.s, params))
}
zkrig <- unlist(krigparts)

thresh <- quantile(zkrig, qthresh)
newfield=as.data.frame(cbind(locs_new_all,zkrig))
colnames(newfield) = c('x','y','z','GRF')
newfield[,'GRF_thresh'] <- ifelse(newfield[,'GRF']>thresh,1,0)

return(newfield)

}
