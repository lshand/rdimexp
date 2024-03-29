% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_exp.R
\name{gen_exp}
\alias{gen_exp}
\title{Generate random image with stationary exponential correlation structure}
\usage{
gen_exp(locs_obs, locs_new, params, qthresh, seed = NA)
}
\arguments{
\item{locs_obs}{- @locs_obs - coordinates of observed data, n x 2 matrix}

\item{locs_new}{- new indices at which to generate new image, N x 2 matrix}

\item{params}{- vector of exponential params, c(sigma2, phi)}

\item{qthresh}{- quantile to threshold at}

\item{seed}{- random seed for generation, default NA}
}
\value{
- prints GRF before and after thresholding (N x 4 matrix)
}
\description{
generates random image of size N x N using parameters fit to observed data of size N x N
}
