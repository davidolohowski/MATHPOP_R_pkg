% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_MATHPOP.R
\name{fit_MATHPOP}
\alias{fit_MATHPOP}
\title{Function that fits the MATHPOP model}
\usage{
fit_MATHPOP(
  Data,
  spat_dom,
  fixed_Theta,
  prior,
  p = 1,
  cf_error,
  M,
  Theta = NULL,
  tune,
  n = 1000,
  prob_model = TRUE,
  seed = 12345,
  burnin = 0.1
)
}
\arguments{
\item{Data}{A data.frame: Each row is an observed GC or point source in the data. Requires at least the columns named \code{x}, \code{y} for the spatial coordinates of GCs (in physical coordinates), and the magnitudes named \code{M}. If there are more columns, they need to be the probabilities a point source is a GC.}

\item{spat_dom}{A List. A list object containing a list called \code{vertices} that gives the vertices of the spatial domain, and a list of the number of integration grid called \code{n_grid}.}

\item{fixed_Theta}{A List. A list that specifies the known parameters of GC system, which contains two list objects \code{gal} and \code{UDG} that specify the respective known parameters of normal galaxies and UDGs. If there are no normal galaxies in the data, \code{gal} does not need to be specified.}

\item{prior}{A List. A list that specifies the parameter values of the prior distributions of the model parameters, specified in a similar fashion as \code{fixed_Theta}.}

\item{p}{A numeric value or vector. Crowding effect. Either a numeric value between \eqn{(0,1)}, or a numeric vector whose entries are all in \eqn{(0,1)} and length equal to \code{n_grid} in \code{spat_dom} with each numeric element being the crowding effect at the location of the spatial grid. In the current implementation, it is default to 1 (no crowding).}

\item{cf_error}{A List. List of parameters for completeness fraction and measurement uncertainties.}

\item{M}{An integer. Total number of iteration to run the MCMC algorithm.}

\item{Theta}{A List. Starting values of the MCMC chain. Default to \code{NULL}, and specified internally.}

\item{tune}{A List. Tuning parameters for initial MCMC pilot run.}

\item{n}{An integer. Initial MCMC pilot run iteration. Default to 1000.}

\item{prob_model}{Logical. Whether the GC data used is a probabilistic catalog or a binary catalog. Default to TRUE.}

\item{seed}{An integer. Random seed value. Default to 12345.}

\item{burnin}{A numeric value. A real number between 0 and 1. Percentage of the sample to be discarded as burn-in for MCMC. Default to 0.1}
}
\value{
A data frame with \code{floor(burn_in*M)} number of rows that gives the posterior sample of the fitted MATHPOP model.
}
\description{
Function that fits the MATHPOP model
}
\keyword{MATHPOP}
\keyword{MCMC}
\keyword{adaptive}
