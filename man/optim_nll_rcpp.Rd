% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{optim_nll_rcpp}
\alias{optim_nll_rcpp}
\title{optim function to estimate phi using profile neg log lik}
\usage{
optim_nll_rcpp(interval, data, ds, prange)
}
\arguments{
\item{interval}{vector of lower and upper bounds for phi}

\item{data}{n x m matrix}

\item{ds}{n x n matrix}

\item{prange}{practical range, if unsure set to max(ds)}
}
\value{
vector with estimated phi and value of obj function
}
\description{
optim function to estimate phi using profile neg log lik
}
