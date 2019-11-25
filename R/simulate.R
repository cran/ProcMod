#' Simulate n points of dimension p.
#'
#' Points are simulated by drawing values of each dimension from a normal
#' distribution of mean 0 and standard deviation equals to 1.
#' The mean of each dimension is forced to 0 (data are centred).
#' By default variable are also scaled to enforce a strandard deviation
#' strictly equal to 1. Covariances between dimensions are not controled.
#' Therefore they are expected to be equal to 0 and reflect only the
#' random distribution of the covariance between two random vectors.
#'
#' @param n an \code{int} value indicating the number of observations.
#' @param p an \code{int} value indicating the number of dimensions (variables)
#'        simulated
#' @param equal_var a \code{logical} value indicating if the dimensions must be scaled
#'        to force \code{sd=1}. \code{TRUE} by default.
#'
#' @return a numeric matrix of \code{n} rows and \code{p} columns
#'
#' @examples
#' sim1 <- simulate_matrix(25,10)
#' class(sim1)
#' dim(sim1)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
simulate_matrix <- function(n, p, equal_var = TRUE) {
  new <- rnorm(n * p, mean = 0, sd = 1)
  dim(new) <- c(n, p)

  new <- scale(new, center = TRUE, scale = equal_var)

  attributes(new)$`scaled:center` <- NULL
  attributes(new)$`scaled:scale` <- NULL
  new.sd <- sqrt(sum(new^2) / (n - 1))
  new / new.sd
}

#' Simulate n points of dimension p correlated to a reference matrix.
#'
#' Simulates a set of point correlated to another set according to the
#' procrustean correlation definition.
#' Points are simulated by drawing values of each dimension from a normal
#' distribution of mean 0 and standard deviation equals to 1.
#' The mean of each dimension is forced to 0 (data are centred).
#' By default variable are also scaled to enforce a strandard deviation
#' strictly equal to 1. Covariances between dimensions are not controled.
#' Therefore they are expected to be equal to 0 and reflect only the
#' random distribution of the covariance between two random vectors.
#' The intensity of the correlation is determined by the \code{r2}
#' parameter.
#'
#' @param reference a numeric matrix to which the simulated data will be correlated
#' @param p an \code{int} value indicating the number of dimensions (variables)
#'        simulated
#' @param r2 the fraction of variation shared between the \code{reference} and the
#'        simulated data
#' @param equal_var a \code{logical} value indicating if the dimensions must be scaled
#'        to force \code{sd=1}. \code{TRUE} by default.
#'
#' @return a numeric matrix of \code{nrow(reference)} rows and \code{p} columns
#'
#' @examples
#' sim1 <- simulate_matrix(25,10)
#' class(sim1)
#' dim(sim1)
#' sim2 <- simulate_correlation(sim1,20,0.8)
#' corls(sim1, sim2)^2
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export

simulate_correlation <- function(reference, p, r2, equal_var = TRUE) {
  n <- nrow(reference)
  maxdim <- max(ncol(reference), p)

  noise <- simulate_matrix(n, p, equal_var = equal_var)

  if (maxdim == p && maxdim > ncol(reference)) {
    # noise is the largest matrix
    YX <- crossprod(noise, reference)
    svd.YX <- svd(YX)
    rot  <- svd.YX$v %*% t(svd.YX$u)
    rotr <- svd.YX$u %*% t(svd.YX$v)

    new = ((reference %*% rot) * sqrt(r2) +
            noise              * sqrt(1 - r2)) %*% rotr
  }
  else {
    # reference is the largest matrix
    YX <- crossprod(reference, noise)
    svd.YX <- svd(YX)
    rot  <- svd.YX$v %*% t(svd.YX$u)
    rotr <- svd.YX$u %*% t(svd.YX$v)

    new = (reference       * sqrt(r2) +
           (noise %*% rot) * sqrt(1 - r2)) %*% rotr
  }

  new <- scale(new, scale = FALSE)
  attributes(new)$`scaled:center` <- NULL
  new.sd <- sqrt(sum(new^2) / (n - 1))
  new / new.sd
}

#simulate_matrix_tree
