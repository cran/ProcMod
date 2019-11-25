#' @include procmod_frame.R
#' @import MASS
#'
NULL

#' Converts a \code{dist} object to a \code{data.frame} object.
#'
#' The created \code{data.frame} has a attribute \code{is.dist} set to
#' the logical value \code{TRUE}.
#'
#' @param x the \code{dist} object to be converted
#' @param row.names NULL or a \code{character} vector giving the row
#'                  names for the data frame. Missing values
#'                  are not allowed.
#' @param optional  logical. If \code{TRUE}, setting row names and converting
#'                  column names (to syntactic names: see make.names)
#'                  is optional. Note that all of R's base package
#'                  as.data.frame() methods use optional only for column
#'                  names treatment, basically with the meaning of
#'                  data.frame(*, check.names = !optional).
#'                  See also the make.names argument of the
#'                  \code{\link[base]{matrix}} method.
#' @param ...       additional arguments to be passed to or from methods.
#'
#' @examples
#' data(bacteria)
#' bacteria_rel_freq <- sweep(bacteria,
#'                            1,
#'                            rowSums(bacteria),
#'                            "/")
#' bacteria_hellinger <- sqrt(bacteria_rel_freq)
#' bacteria_dist <- dist(bacteria_hellinger)
#' bdf <- as.data.frame(bacteria_dist)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
as.data.frame.dist <- function(x, row.names = NULL, optional = FALSE, ...) {

  x <- as.data.frame(
    as.matrix(x),
    row.names = row.names,
    optional = optional,
    ...)

  attr(x, "is.dist") <- TRUE

  x
}


#' Project a distance matrix in a euclidean space (NMDS).
#'
#' Project a set of points defined by a distance matrix in
#' an eucleadean space using the Kruskal's Non-metric
#' Multidimensional Scaling. This function is mainly a simplified
#' interface on the \code{\link[MASS]{isoMDS}} function using as
#' much as possible dimensions to limit the stress. The aims of this
#' NDMS being only to project point in an orthogonal space therefore
#' without any correlation between axis. Because a non-metric method
#' is used no condition is required on the used distance.
#'
#' @param distances   a \code{\link[stats]{dist}} object or a
#'                    \code{\link[base]{matrix}}
#'                    object representing a distance matrix.
#' @param maxit       The maximum number of iterations.
#' @param trace       Logical for tracing optimization. Default \code{TRUE}.
#' @param tol         convergence tolerance.
#' @param p           Power for Minkowski distance in the configuration space.
#'
#' @return a numeric matrix with at most \code{n-1} dimensions, with
#'         \code{n} the number pf observations. This matrix defines the
#'         coordinates of each point in the orthogonal space.
#'
#' @examples
#' data(bacteria)
#' bacteria_rel_freq <- sweep(bacteria,
#'                            1,
#'                            rowSums(bacteria),
#'                            "/")
#' bacteria_hellinger <- sqrt(bacteria_rel_freq)
#' bacteria_dist <- dist(bacteria_hellinger)
#'
#' project <- nmds(bacteria_dist)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
nmds <- function(distances,
                 maxit = 100, trace = FALSE,
                 tol = 0.001, p = 2) {


  if (inherits(distances, "matrix")) {
    distances <- as.dist(distances)
  }

  stopifnot(inherits(distances, "dist"))

  k <- attributes(distances)$Size - 1

  y <- suppressWarnings(cmdscale(distances, k))
  k <- ncol(y)

  n <- isoMDS(distances,
    y = y,
    k = k,
    maxit = maxit,
    trace = trace,
    tol = tol,
    p = p
  )

  p <- n$points
  attributes(p)$stress <- n$stress
  attributes(p)$projected <- "nmds"

  p
}

#' Project a distance matrix in a euclidean space (PCOA).
#'
#' Project a set of points defined by a distance matrix in
#' an eucleadean space using the Principal Coordinates Analysis
#' method. This function is mainly a simplified
#' interface on the \code{\link[stats]{cmdscale}} function using as
#' much as possible dimensions for the projection. The aims of this
#' PCoA being only to project point in an orthogonal space therefore
#' without any correlation between axis. Because a metric method
#' is used the used distance must be euclidean
#' (cf \code{\link[ProcMod]{is_euclid}}).
#'
#' @param distances   a \code{\link[stats]{dist}} object or a
#'                    \code{\link[base]{matrix}}
#'                    object representing a distance matrix.
#'
#' @return a numeric matrix with at most \code{n-1} dimensions, with
#'         \code{n} the number pf observations. This matrix defines the
#'         coordinates of each point in the orthogonal space.
#'
#' @examples
#' data(bacteria)
#' bacteria_rel_freq <- sweep(bacteria,
#'                            1,
#'                            rowSums(bacteria),
#'                            "/")
#' bacteria_hellinger <- sqrt(bacteria_rel_freq)
#' bacteria_dist <- dist(bacteria_hellinger)
#'
#' project <- pcoa(bacteria_dist)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
pcoa <- function(distances) {

  if (inherits(distances, "matrix")) {
    distances <- as.dist(distances)
  }

  stopifnot(inherits(distances, "dist"))

  k <- attributes(distances)$Size - 1
  x <- suppressWarnings(cmdscale(distances, k = k))
  attributes(x)$projected <- "pcoa"

  x
}

#' Project a set of points in a euclidean space (PCA).
#'
#' Project a set of points defined by a set of numeric variables in
#' an eucleadean space using the pricipal componant analysis.
#' This function is mainly a simplified
#' interface on the \code{\link[stats]{prcomp}} function using as
#' much as possible dimensions to keep all the variation. The aims of this
#' PCA being only to project point in an orthogonal space therefore
#' without any correlation between axis. Data are centered by not scaled by
#' default.
#'
#' @param data a numeric matrix describing the points
#' @param scale a \code{logical} value indicating if the dimensions must be scaled
#'        to force for every column that \code{sd=1}. \code{FALSE} by default.
#'
#' @return a numeric matrix with at most \code{n-1} dimensions, with
#'         \code{n} the number pf observations. This matrix defines the
#'         coordinates of each point in the orthogonal space.
#'
#' @examples
#' data(bacteria)
#' bacteria_rel_freq <- sweep(bacteria,
#'                            1,
#'                            rowSums(bacteria),
#'                            "/")
#' bacteria_hellinger <- sqrt(bacteria_rel_freq)
#'
#' project <- pca(bacteria_hellinger)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
pca <- function(data, scale = FALSE) {

  p <- prcomp(data,
    retx = TRUE,
    center = TRUE,
    scale. = scale,
    tol = 0
  )

  x <- p$x
  attributes(x)$projected <- "pca"

  x
}

#' Double centering of a matrix.
#'
#' colSums and rowSums of the returned matrix are all equal to zero.
#'
#' Inspired from the algorithm described in stackoverflow
#' \url{https://stackoverflow.com/questions/43639063/double-centering-in-r}
#'
#' @param m a \code{numeric} matrix
#'
#' @return a \code{numeric} matrix centred by rows
#'           and columns
#'
#' @examples
#' data(bacteria)
#' bact_bc <- bicenter(bacteria)
#' sum(rowSums(bact_bc))
#' sum(colSums(bact_bc))
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
bicenter <- function(m) {

  # compute the row-wise and column-wise mean matrices
  m0 <- m * 0
  R <- m0 + rowMeans(m)
  C <- t(m0 + colMeans(m))

  # substract them and add the grand mean
  m - R - C + mean(m[])
}

#' Test if the distance matrix is euclidean.
#'
#' Actually a simplified version of the ADE4 implementation
#' (\code{\link[ade4]{is.euclid}}).
#'
#' @param distances	an object of class 'dist'
#' @param tol	a tolerance threshold : an eigenvalue is
#'        considered positive if it is larger than
#'        -tol*lambda1 where lambda1 is the largest eigenvalue.
#'
#' @examples
#' library(vegan)
#' data(bacteria)
#'
#' bacteria_rel_freq <- sweep(bacteria,
#'                            1,
#'                            rowSums(bacteria),
#'                            "/")
#'
#' bacteria_bray <- vegdist(bacteria_rel_freq,method = "bray")
#' is_euclid(bacteria_bray)
#'
#' bacteria_chao <- vegdist(floor(bacteria*10000),method = "chao")
#' is_euclid(bacteria_chao)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
is_euclid <- function(distances, tol = 1e-07) {
  if (!inherits(distances, "dist")) {
    stop("Object of class 'dist' expected")
  }

  if (any(distances < tol)) {
    warning("Zero distance(s)")
  }

  distances <- as.matrix(distances)
  n <- ncol(distances)
  delta <- -0.5 * ProcMod::bicenter(distances * distances)
  lambda <- eigen(delta, symmetric = TRUE, only.values = TRUE)$values
  w0 <- lambda[n] / lambda[1]

  w0 > -tol
}

is_orthogonal <- function(x) {
  stopifnot(is_procmod_frame(x))

  !is.null(attr(x, "projected"))
}

#' Project a dataset in a euclidean space.
#'
#' Project a set of points defined by a distance matrix
#' or a set of variables in an eucleadean space.
#' If the distance matrix is a metric, this is done using
#' the \code{\link[ProcMod]{pcoa}} function,
#' for other distance the \code{\link[ProcMod]{nmds}} is used.
#' When points are described by a set of variable the
#' \code{\link[ProcMod]{nmds}} is used.
#'
#' @param data a numeric matrix describing the points
#' @param tol	a tolerance threshold : an eigenvalue is
#'        considered positive if it is larger than
#'        -tol*lambda1 where lambda1 is the largest eigenvalue.
#' @param scale a \code{logical} value indicating if the dimensions must be scaled
#'        to force for every column that \code{sd=1}. \code{FALSE} by default.
#'
#' @param ... other parameters specific to some implementation
#'
#' @return a numeric matrix with at most \code{n-1} dimensions, with
#'         \code{n} the number pf observations. This matrix defines the
#'         coordinates of each point in the orthogonal space.
#'
#' @examples
#' library(vegan)
#' data(bacteria)
#' data(eukaryotes)
#' data(soil)
#'
#' dataset <- procmod_frame(euk = vegdist(decostand(eukaryotes,
#'                                                  method = "hellinger"),
#'                                        method = "euclidean"),
#'                          bac = vegdist(decostand(bacteria,
#'                                                  method = "hellinger"),
#'                                        method = "euclidean"),
#'                          soil = scale(soil,
#'                                       center = TRUE,
#'                                       scale  = TRUE))
#'
#' dp <- ortho(dataset)
#'
#' bacteria_rel_freq <- sweep(bacteria,
#'                            1,
#'                            rowSums(bacteria),
#'                            "/")
#' bacteria_hellinger <- sqrt(bacteria_rel_freq)
#' bacteria_dist <- dist(bacteria_hellinger)
#'
#' project <- ortho(bacteria_dist)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
ortho <- function(data, ...) {
  UseMethod("ortho", data)
}

#' @rdname ortho
#' @export
ortho.dist <- function(data, tol = 1e-7, ...) {
  if (ProcMod::is_euclid(data, tol = tol)) {
    return(ProcMod::pcoa(data))
  }

  ProcMod::nmds(data)
}

#' @rdname ortho
#' @export
ortho.matrix <- function(data, scale = FALSE, ...) {
  pca(data, scale = scale)
}

#' @rdname ortho
#' @export
ortho.data.frame <- function(data, scale = FALSE, ...) {
  if (!is.null(attributes(data)$is.dist) &&
    attributes(data)$is.dist == TRUE) {
    return(ortho(as.dist(as.matrix(data))))
  }

  pca(as.matrix(data), scale = scale)
}


#' @rdname ortho
#' @export
ortho.procmod_frame <- function(data, ...) {
  if (is_orthogonal(data)) {
    return(data)
  }

  n <- ncol(data)

  p <- vector(mode = "character", length = n)
  for (i in seq_len(n)) {
    xt <- ortho(data[[i]])
    p[i] <- attributes(xt)$projected
    data[[i]] <- xt
  }

  names(p) <- names(data)
  attributes(data)$projected <- p

  data
}
