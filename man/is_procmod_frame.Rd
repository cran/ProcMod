% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/procmod_frame.R
\name{is_procmod_frame}
\alias{is_procmod_frame}
\title{Check if an object is a ProcMod Frame.}
\usage{
is_procmod_frame(x)
}
\arguments{
\item{x}{a R \code{object to test}}
}
\value{
a \code{logical} value equals to
        \code{TRUE} if \code{x} is a \code{procmod_frame},
        \code{FALSE} otherwise.
}
\description{
Check if an object is a ProcMod Frame.
}
\examples{
# Builds a procmod_frame with two random matrices
m1 <- simulate_matrix(10,20)
m2 <- simulate_matrix(10,30)
pmf <- procmod_frame(m1 = m1, m2 = m2)

# Returns TRUE
is_procmod_frame(pmf)

# Returns FALSE
is_procmod_frame(3)

}
\author{
Eric Coissac

Christelle Gonindard-Melodelima
}
