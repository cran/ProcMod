#' @include internals.R
#'
NULL

#' Rotate the \code{src} matrix to fit into the space of the \code{dest} matrix.
#'
#' The optimal rotation is computed according to the procruste methode.
#' Rotation is based on singular value decomposition (SVD).
#' No scaling and no centrering are done, before computing the SVD.
#'
#' @param src a numeric matrix to be rotated
#' @param dest a numeric matrix used as reference space
#'
#' @return a numeric matrix
#'
#' @examples
#' # Generates two random matrices of size 10 x 15
#' m1 <- simulate_matrix(10, 15)
#' m2 <- simulate_matrix(10, 20)
#'
#' # Rotates matrix m1 on m2
#' mr <- protate(m1, m2)
#'
#' @author Christelle Gonindard-Melodelima
#' @author Eric Coissac
#' @export
protate <- function(src, dest) {
  YX <- crossprod(dest, src)
  svd.YX <- svd(YX)
  rot <- svd.YX$v %*% t(svd.YX$u)
  src %*% rot
}
