#' compute autocorrelation matrix of img using fast fourier transforms
#'
#'
#' @param x - n x n raw image matrix
#' @return - n x n autocorrelation matrix
#' @export

autocorr <- function(x){
  fftx = fft(x)
  fftx2 = Re(signal::ifft(fftx*Conj(fftx)))/prod(dim(fftx))

  shift <- function(x){
    nc = dim(x)[2]
    order = 1:nc
    neworder = sapply(order,function(x) x+as.integer(nc/2+1))
    neworder = ifelse(neworder > nc, neworder-nc,neworder)
    return(x[,neworder])
  }

  fftx3 = shift(apply(fftx2,2,fftshift))
  return(fftx3)
}

#' compute radial profile of raw image
#'
#'
#' @param locs x 2 matrix of locations [x,y] #column name matters
#' @param x x 1 vector of observed value
#' @return - list of 2: vector of distances and vector of respective radial profile values
#' @export

radial_profile <-function(locs, x, sigd=0){ #returns correlation matrix

  x2=autocorr(matrix(x,nrow=length(unique(locs[,'x'])),ncol=length(unique(locs[,'y']))))
  center = apply(locs,2,mean)
  inds1 = locs[,1]
  inds2 = locs[,2]
  r = sqrt((inds1-center[1])^2+(inds2-center[2])^2) #euclidean distance
  center_value = x2[dim(x2)[1]*0.5,dim(x2)[2]*0.5]

  if(sigd==0){
    rfac = as.factor(as.integer(r))
  }else{
    rfac = as.factor(round(r,sigd))
  }

  tbin = freqweight(rfac,as.vector(round(x2,4)))
  nr = as.vector(table(rfac))

  rprof=tbin/nr
  rpmax <- which.max(rprof)
  if(rpmax==length(rprof)){
    rprof <- c(rprof,center_value)
    rvec=sort(c(0,unique(as.numeric(rfac))),decreasing = TRUE)
  }else{
    rvec=sort(c(0,unique(as.numeric(rfac))),decreasing = FALSE)
    rprof=c(center_value,tbin/nr)
  }

  return(list(rvec=rvec,covvec=rprof))
}
