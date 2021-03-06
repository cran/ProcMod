% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covls.R
\name{.Trace}
\alias{.Trace}
\title{Compute the trace of a square matrix.}
\usage{
.Trace(X)
}
\arguments{
\item{X}{a square matrix}
}
\value{
the trace of X
}
\description{
The trace of a square matrix is defined as the sum
of its diagonal elements.
}
\note{
Internal function do not use.
}
\examples{
m <- matrix(1:16, nrow = 4)
ProcMod:::.Trace(m)
}
\author{
Eric Coissac

Christelle Gonindard-Melodelima
}
