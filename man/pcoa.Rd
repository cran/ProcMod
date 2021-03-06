% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate.R
\name{pcoa}
\alias{pcoa}
\title{Project a distance matrix in a euclidean space (PCOA).}
\usage{
pcoa(distances)
}
\arguments{
\item{distances}{a \code{\link[stats]{dist}} object or a
\code{\link[base]{matrix}}
object representing a distance matrix.}
}
\value{
a numeric matrix with at most \code{n-1} dimensions, with
        \code{n} the number pf observations. This matrix defines the
        coordinates of each point in the orthogonal space.
}
\description{
Project a set of points defined by a distance matrix in
an eucleadean space using the Principal Coordinates Analysis
method. This function is mainly a simplified
interface on the \code{\link[stats]{cmdscale}} function using as
much as possible dimensions for the projection. The aims of this
PCoA being only to project point in an orthogonal space therefore
without any correlation between axis. Because a metric method
is used the used distance must be euclidean
(cf \code{\link[ProcMod]{is_euclid}}).
}
\examples{
data(bacteria)
bacteria_rel_freq <- sweep(bacteria,
                           1,
                           rowSums(bacteria),
                           "/")
bacteria_hellinger <- sqrt(bacteria_rel_freq)
bacteria_dist <- dist(bacteria_hellinger)

project <- pcoa(bacteria_dist)

}
\author{
Eric Coissac

Christelle Gonindard-Melodelima
}
