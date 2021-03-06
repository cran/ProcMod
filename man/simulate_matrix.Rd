% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate_matrix}
\alias{simulate_matrix}
\title{Simulate n points of dimension p.}
\usage{
simulate_matrix(n, p, equal_var = TRUE)
}
\arguments{
\item{n}{an \code{int} value indicating the number of observations.}

\item{p}{an \code{int} value indicating the number of dimensions (variables)
simulated}

\item{equal_var}{a \code{logical} value indicating if the dimensions must be scaled
to force \code{sd=1}. \code{TRUE} by default.}
}
\value{
a numeric matrix of \code{n} rows and \code{p} columns
}
\description{
Points are simulated by drawing values of each dimension from a normal
distribution of mean 0 and standard deviation equals to 1.
The mean of each dimension is forced to 0 (data are centred).
By default variable are also scaled to enforce a strandard deviation
strictly equal to 1. Covariances between dimensions are not controled.
Therefore they are expected to be equal to 0 and reflect only the
random distribution of the covariance between two random vectors.
}
\examples{
sim1 <- simulate_matrix(25,10)
class(sim1)
dim(sim1)

}
\author{
Eric Coissac

Christelle Gonindard-Melodelima
}
