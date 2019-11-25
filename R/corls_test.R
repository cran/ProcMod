#' @include covls.R
#' @import permute
#' @author Christelle Gonindard-Melodelima
#' @author Eric Coissac
NULL

#' Generate permutation matrix according to a schema.
#'
#' The permutation schema is defined using the `how` function.
#' The implementation of this function is inspired
#' from the VEGAN package and reproduced here to avoid an extra
#' dependency on an hidden vegan function.
#'
#' @param permutations a list of control values for the permutations as returned
#'              by the function \code{\link[permute]{how}}, or the number of
#'              permutations required.
#' @param n numeric; the number of observations in the sample set.
#'              May also be any object that nobs knows about;
#'              see \code{\link[permute]{nobs}} methods.
#' @param strata A factor, or an object that can be coerced to a
#'               factor via as.factor, specifying the strata for permutation.
#'
#' @note Internal function do not use.
#'
#' @rdname internal.getPermuteMatrix
.getPermuteMatrix = function(permutations, n, strata = NULL)
{
  if (length(permutations) == 1) {
    permutations <- permute::how(nperm = permutations)
  }
  if (!missing(strata) && !is.null(strata)) {
    if (inherits(permutations, "how") && is.null(permute::getBlocks(permutations)))
      setBlocks(permutations) <- strata
  }
  if (inherits(permutations, "how"))
    permutations <- permute::shuffleSet(n, control = permutations)
  else {
    if (!is.integer(permutations) && !all(permutations == round(permutations)))
      stop("permutation matrix must be strictly integers: use round()")
  }
  if (is.null(attr(permutations, "control")))
    attr(permutations, "control") <- structure(list(within = list(type = "supplied matrix"),
                                            nperm = nrow(permutations)), class = "how")
  permutations
}



#' Monte-Carlo Test on the sum of the singular values of a procustean rotation.
#'
#' performs a Monte-Carlo Test on the sum of the singular values of a
#' procustean rotation (see \code{\link[ade4]{procuste.rtest}}).
#'
#' @param ...   the set of matrices or a \code{\link[ProcMod]{procmod_frame}}
#'              object.
#' @param permutations a list of control values for the permutations as returned
#'              by the function \code{\link[permute]{how}}, or the number of
#'              permutations required.
#' @param p_adjust_method the multiple test correction method used
#'              to adjust p values. \code{p_adjust_method}
#'              belongs one of the folowing values: \code{"holm"},
#'              \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'              \code{"BH"}, \code{"BY"}, \code{"fdr"},  \code{"none"}.
#'              The default is,set to \code{"holm"}.
#'
#' @references {
#'  \insertRef{Jackson:95:00}{ProcMod}
#' }
#'
#' @examples
#' A <- simulate_matrix(10,3)
#' B <- simulate_matrix(10,5)
#' C <- simulate_correlation(B,10,r2=0.6)
#'
#' # Computes the correlation matrix
#' data <- procmod_frame(A = A, B = B, C = C)
#'
#' corls_test(A, B, permutations = 100)
#' corls_test(B, C, permutations = 100)
#' corls_test(data, permutations = 100)
#'
#' @seealso \code{\link[stats]{p.adjust}}
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
corls_test <- function(...,
                       permutations = permute::how(nperm = 999),
                       p_adjust_method="holm") {
  eps <- sqrt(sqrt(.Machine$double.eps))
  xs <- list(...)

  if (length(xs) == 1) {
    x <- xs[[1]]
    if (is_procmod_frame(x)) {
      xs <- x
    } else {
      xs <- procmod_frame(x)
    }
  }
  else {
    xs <- as_procmod_frame(xs)
  }

  x_names <- names(xs)

  xs <- ortho(xs)


  cov <- varls(xs, nrand = 0)
  lcov <- cov - eps
  ngreater <- array(0,dim = dim(cov))

  n <- nrow(xs)
  nx <- length(xs)
  pmatrix <- .getPermuteMatrix(permutations, n)

  if (ncol(pmatrix) != n) {
    stop(gettextf(
      "'permutations' have %d columns, but data have %d observations",
      ncol(pmatrix), n
    ))
  }

  npermutation <- nrow(pmatrix)
  for (i in seq_len(npermutation)) {
    ps <- sample(1:npermutation,
      size = nx,
      replace = FALSE
    )

    rcov = varls(as_procmod_frame(
        lapply(
          1:nx,
          function(j) xs[[j]][pmatrix[ps[j], ], ]
        )
      ),
      nrand = 0
      )
    ngreater <- ngreater + (
       rcov >= lcov)
  }

  p_values <- ngreater / npermutation
  diag(p_values) <- 0

  c_p_values <- p.adjust(p_values[upper.tri(p_values,diag = FALSE)],
           method = p_adjust_method,
           n = (nx - 1) * nx / 2)

  p_values[upper.tri(p_values,diag = FALSE)] <- c_p_values
  p_values <- as.matrix(Matrix::forceSymmetric(p_values, uplo = "U"))
  colnames(p_values) <- x_names
  rownames(p_values) <- x_names

  p_values

}
