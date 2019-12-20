#' Fit Dimension Expansion model to data and regenerate random image of same size
#'
#'
#' @param locs_obs - coordinates of observed data, N x 2 matrix
#' @param datamat - N x m matrix of m replicates of size n at observed locations
#' @param per - proportion of points n to uniformly sample across image, default p=0.05
#' @param opt_param_init - initial n+2 parameter estimates for theta. If default NA specified, MLE and rnorm(0,1) is used.
#' @param lambdas - list specify tuning parameters lambda1 and lambda2. If NA, grid search performed over grid_lambdas, default: NA
#' @param grid_lambdas - grid of lambdas to optimize over, list of 2 vectors. default: list(lambda1=c(0.01,0.1,1,10),lambda2=c(0.0001,0.001,0.01,0.1))
#' @param seed - default NA
#' @param extra - if TRUE prints extra info such as mlconv matrix
#' @param progress - if TRUE (default) print progress bar
#' @return - list including estimated theta, optimal lambdas and @locs_obs
#' @export

fit_dimex<-function(locs_obs, datamat, per=0.05, opt_param_init=NA, lambdas=list(lambda1=NA, lambda2=NA), grid_lambdas=list(lambda1=c(0.01,0.1,1,10),lambda2=c(0.0001,0.001,0.01,0.1)), seed=NA, extra=FALSE, progress=TRUE){

  if(!is.na(lambdas$lambda1)){
  maxit = 2
  if(progress) pb <- txtProgressBar(min = 0, max = maxit, style = 3)
  }else{
  nl1 <- length(grid_lambdas$lambda1)
  nl2 <- length(grid_lambdas$lambda2)
  maxit = nl1 * nl2 + 1
  if(progress) pb <- txtProgressBar(min = 0, max = maxit, style = 3)
  }

  if(is.na(seed)){
    seed = sample(1:1000)
  }
  set.seed(seed)

  ns=ceiling(dim(datamat)[1]*per)
  obs <- sample(1:dim(datamat)[1],ns)
  locs_obs_clip <- as.matrix(locs_obs[obs,1:2])
  datamat_clip<-datamat[obs,]
  ds_clip <- dist_euclid(locs_obs_clip)
  CovObs_clip <- cov(t(datamat_clip)) #construct covariance matrix C(ds_clip)
  empvario_clip <- emp_vario(CovObs_clip) #variogram

  if(length(opt_param_init) != (ns+2)){
    param_init <- optim_nll_rcpp(interval=c(0,200),data=datamat_clip,ds=ds_clip,prange=max(ds_clip))[1] #decay parameter phi
    sig_init <- mean(diag(t(datamat_clip)%*%solve(exp(-param_init*ds_clip),datamat_clip))/ns)
    locs_int = cbind(locs_obs_clip,rnorm(ns,0,mean(apply(locs_obs_clip,2,sd))))
    opt_param_init = c(sig_init,param_init,locs_int[,3])
	
	if(progress){
    Sys.sleep(0.1)
    setTxtProgressBar(pb, 1)# update progress bar
	}
  }else{
    locs_int = cbind(locs_obs_clip,opt_param_init[-(1:2)])
	
	if(progress){
    Sys.sleep(0.1)
    setTxtProgressBar(pb, 1)# update progress bar
	}
  }


  fit_dimex1 <-function(lambda1){
    opt_params_tmp = optim_z_rcpp(c(log(opt_param_init[1:2]),opt_param_init[-(1:2)]),locs_obs=locs_obs_clip,emp_vario=empvario_clip, lambda=lambda1)
    #opt_params = optim_results$par
    opt_params_tmp[1:2] <- exp(opt_params_tmp[1:2])
    return(opt_params_tmp)
    }


  if(!is.na(lambdas$lambda1)){
    opt_params <- fit_dimex1(lambdas$lambda1)
    opt_lambdas <- lambdas
	
	if(progress){
    Sys.sleep(0.1)
    setTxtProgressBar(pb, 2)# update progress bar
	}
  }else{
    lpairs <- expand.grid(grid_lambdas$lambda2, grid_lambdas$lambda1)
    mlconv <- cbind(lpairs[,c(2,1)],rep(0,dim(lpairs)[1]))

    locs_int1 <- lapply(grid_lambdas$lambda1,function(x) fit_dimex1(x))

    for(i in 1:nl1){
      locs_int[,3]=locs_int1[[i]][-(1:2)]

      for(j in 1:nl2){
        row = (i-1)*nl2+j
        fit.s = fields::Tps(locs_int[,-3], locs_int[,3],lambda=grid_lambdas$lambda2[j],give.warnings=FALSE)
        locs_int[,3] = predict(fit.s,locs_int[,-3], locs=0)
        ds=dist_euclid(as.matrix(locs_int))
        #ml_est_lrnd<-optimize(prof_nll,c(0.01,100),data=datamat, ds=ds, prange=max(ds))
        ml_est_lrnd<-optim_nll_rcpp(interval=c(0.01,100),data=datamat_clip, ds=ds, prange=max(ds))
        mlconv[row,3]<- ml_est_lrnd[2]
		
		if(progress){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, row+1)# update progress bar
		}
      }
    }

    opt_lambdas <- as.numeric(mlconv[min(which.min(mlconv[,3])),1:2])
    opt_params <- locs_int1[[which(grid_lambdas$lambda1==opt_lambdas[1])]]
  }

  if(extra){
    return(list(opt_params=opt_params, lambdas=opt_lambdas, locs_obs=locs_obs_clip, mlconv=mlconv))
  }else{
    return(list(opt_params=opt_params, lambdas=list(lambda1 <- opt_lambdas[1], lambda2 <- opt_lambdas[2]), locs_obs=locs_obs_clip))
  }

}
