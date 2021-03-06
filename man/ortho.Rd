% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate.R
\name{ortho}
\alias{ortho}
\alias{ortho.dist}
\alias{ortho.matrix}
\alias{ortho.data.frame}
\alias{ortho.procmod_frame}
\title{Project a dataset in a euclidean space.}
\usage{
ortho(data, ...)

\method{ortho}{dist}(data, tol = 1e-07, ...)

\method{ortho}{matrix}(data, scale = FALSE, ...)

\method{ortho}{data.frame}(data, scale = FALSE, ...)

\method{ortho}{procmod_frame}(data, ...)
}
\arguments{
\item{data}{a numeric matrix describing the points}

\item{...}{other parameters specific to some implementation}

\item{tol}{a tolerance threshold : an eigenvalue is
considered positive if it is larger than
-tol*lambda1 where lambda1 is the largest eigenvalue.}

\item{scale}{a \code{logical} value indicating if the dimensions must be scaled
to force for every column that \code{sd=1}. \code{FALSE} by default.}
}
\value{
a numeric matrix with at most \code{n-1} dimensions, with
        \code{n} the number pf observations. This matrix defines the
        coordinates of each point in the orthogonal space.
}
\description{
Project a set of points defined by a distance matrix
or a set of variables in an eucleadean space.
If the distance matrix is a metric, this is done using
the \code{\link[ProcMod]{pcoa}} function,
for other distance the \code{\link[ProcMod]{nmds}} is used.
When points are described by a set of variable the
\code{\link[ProcMod]{nmds}} is used.
}
\examples{
library(vegan)
data(bacteria)
data(eukaryotes)
data(soil)

dataset <- procmod_frame(euk = vegdist(decostand(eukaryotes,
                                                 method = "hellinger"),
                                       method = "euclidean"),
                         bac = vegdist(decostand(bacteria,
                                                 method = "hellinger"),
                                       method = "euclidean"),
                         soil = scale(soil,
                                      center = TRUE,
                                      scale  = TRUE))

dp <- ortho(dataset)

bacteria_rel_freq <- sweep(bacteria,
                           1,
                           rowSums(bacteria),
                           "/")
bacteria_hellinger <- sqrt(bacteria_rel_freq)
bacteria_dist <- dist(bacteria_hellinger)

project <- ortho(bacteria_dist)

}
\author{
Eric Coissac

Christelle Gonindard-Melodelima
}
