% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate_correlation}
\alias{simulate_correlation}
\title{Simulate n points of dimension p correlated to a reference matrix.}
\usage{
simulate_correlation(reference, p, r2, equal_var = TRUE)
}
\arguments{
\item{reference}{a numeric matrix to which the simulated data will be correlated}

\item{p}{an \code{int} value indicating the number of dimensions (variables)
simulated}

\item{r2}{the fraction of variation shared between the \code{reference} and the
simulated data}

\item{equal_var}{a \code{logical} value indicating if the dimensions must be scaled
to force \code{sd=1}. \code{TRUE} by default.}
}
\value{
a numeric matrix of \code{nrow(reference)} rows and \code{p} columns
}
\description{
Simulates a set of point correlated to another set according to the
procrustean correlation definition.
Points are simulated by drawing values of each dimension from a normal
distribution of mean 0 and standard deviation equals to 1.
The mean of each dimension is forced to 0 (data are centred).
By default variable are also scaled to enforce a strandard deviation
strictly equal to 1. Covariances between dimensions are not controled.
Therefore they are expected to be equal to 0 and reflect only the
random distribution of the covariance between two random vectors.
The intensity of the correlation is determined by the \code{r2}
parameter.
}
\examples{
sim1 <- simulate_matrix(25,10)
class(sim1)
dim(sim1)
sim2 <- simulate_correlation(sim1,20,0.8)
corls(sim1, sim2)^2

}
\author{
Eric Coissac

Christelle Gonindard-Melodelima
}
