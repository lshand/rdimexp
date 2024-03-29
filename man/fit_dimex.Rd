% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_dimex.R
\name{fit_dimex}
\alias{fit_dimex}
\title{Fit Dimension Expansion model to data and regenerate random image of same size}
\usage{
fit_dimex(locs_obs, datamat, per = 0.05, opt_param_init = NA,
  lambdas = list(lambda1 = NA, lambda2 = NA),
  grid_lambdas = list(lambda1 = c(0.01, 0.1, 1, 10), lambda2 = c(1e-04,
  0.001, 0.01, 0.1)), seed = NA, extra = FALSE, progress = TRUE)
}
\arguments{
\item{locs_obs}{- coordinates of observed data, N x 2 matrix}

\item{datamat}{- N x m matrix of m replicates of size n at observed locations}

\item{per}{- proportion of points n to uniformly sample across image, default p=0.05}

\item{opt_param_init}{- initial n+2 parameter estimates for theta. If default NA specified, MLE and rnorm(0,1) is used.}

\item{lambdas}{- list specify tuning parameters lambda1 and lambda2. If NA, grid search performed over grid_lambdas, default: NA}

\item{grid_lambdas}{- grid of lambdas to optimize over, list of 2 vectors. default: list(lambda1=c(0.01,0.1,1,10),lambda2=c(0.0001,0.001,0.01,0.1))}

\item{seed}{- default NA}

\item{extra}{- if TRUE prints extra info such as mlconv matrix}

\item{progress}{- if TRUE (default) print progress bar}
}
\value{
- list including estimated theta, optimal lambdas and @locs_obs
}
\description{
Fit Dimension Expansion model to data and regenerate random image of same size
}
