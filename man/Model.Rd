% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-class.R
\name{Model}
\alias{Model}
\title{An S4 super class to represent a regression-like model.}
\usage{
Model(
  ndocs,
  nchain = 1,
  mu0 = NaN,
  sigma0 = NaN,
  eta_start = NaN,
  eta = NaN,
  loglike = NaN,
  logpost = NaN,
  waic = NaN,
  se_waic = NaN,
  p_eff = NaN,
  lpd = NaN
)

Model(
  ndocs,
  nchain = 1,
  mu0 = NaN,
  sigma0 = NaN,
  eta_start = NaN,
  eta = NaN,
  loglike = NaN,
  logpost = NaN,
  waic = NaN,
  se_waic = NaN,
  p_eff = NaN,
  lpd = NaN
)
}
\description{
An S4 super class to represent a regression-like model.

Helper function (constructor) for Model class.
}
\section{Slots}{

\describe{
\item{\code{ndocs}}{The number of documents/observations.}

\item{\code{nchain}}{The number of iterations of the Gibbs sampler.}

\item{\code{mu0}}{A (p + 1) x 1 matrix of prior means for eta.}

\item{\code{sigma0}}{A (p + 1) x (p + 1) prior covariance matrix for eta.}

\item{\code{eta_start}}{A (p + 1) x 1 matrix of starting values for eta.}

\item{\code{eta}}{A nchain x (p + 1) matrix of draws of regression coefficients.}

\item{\code{loglike}}{A nchain x 1 vector of the log-likelihood (up to an additive
constant).}

\item{\code{logpost}}{A nchain x 1 vector of the log-posterior (up to an additive
constant).}

\item{\code{waic}}{WAIC (up to an additive constant) on the deviance scale.}

\item{\code{se_waic}}{Standard error of the WAIC.}

\item{\code{p_eff}}{The effective number of parameters.}

\item{\code{lpd}}{A nchain x ndocs matrix of predictive posterior likelihoods (NOT
log-likelihoods).}

\item{\code{extra}}{A list of additional model fitting information. Contains
time_elapsed, start_time, end_time, corrected_label_switching, and call.}
}}

