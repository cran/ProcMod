#' @include procmod.R
#' @include procmod_frame.R
#' @include multivariate.R
#' @import foreach
#' @import stats
#'
#' @author Christelle Gonindard-Melodelima
#' @author Eric Coissac
NULL

.has_doParallel <- is.element("doParallel",installed.packages())
if (.has_doParallel) require(doParallel)

#' Compute the trace of a square matrix.
#'
#' The trace of a square matrix is defined as the sum
#' of its diagonal elements.
#'
#' @param X a square matrix
#' @return the trace of X
#'
#' @examples
#' m <- matrix(1:16, nrow = 4)
#' ProcMod:::.Trace(m)
#' @note Internal function do not use.
#'
#' @rdname internal.Trace
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#'
.Trace <- function(X) sum(diag(X))

#' Convert a covariance matrix to a correlation matric
#'
#' @param c a square covariance/variance matrix
#' @return correlation matrix corresponding to \code{c}
#'
#' @note Internal function do not use.
#'
#' @rdname internal.var2cor
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#'
.var2cor <- function(c) {
  v <- sqrt(diag(c))
  vv <- v %o% v
  c / vv
}

#' Procrustean Correlation, and Variance / Covariance Matrices.
#'
#' \code{varls}, \code{corls}, \code{corls_partial} compute the procrustean
#' variance / covariance, correlation, or partial correlation matrices
#' between a set of real matrices and \code{\link[stats]{dist}} objects.
#'
#' Procrustean covariance between two matrices X and Y, is defined as the sum
#' of the singular values of the X'Y matrix \insertCite{Gower:71:00,Lingoes:74:00}{ProcMod}.
#' Both the X and Y matrices must have the same number of rows.
#'
#' The variances and covariances and correlations are corrected
#' to avoid over fitting \insertCite{Coissac-Eric:19:00}{ProcMod}.
#'
#' Partial correlation coefficients are computed by inverting the correlation followed
#' by a normalisation by the diagonal of the inverted matrix.
#'
#' The inputs must be numeric matrices or \code{\link[stats]{dist}} object.
#' The set of input matrices can be aggregated un a
#' \code{\link[ProcMod]{procmod_frame}}.
#'
#' Before computing the coefficients, matrices are projected into an
#' orthogonal space using the \code{\link[ProcMod]{ortho}} function.
#'
#' The denominator n - 1 is used which gives an unbiased estimator of the
#' (co)variance for i.i.d. observations.
#'
#' Scaling a covariance matrix into a correlation one can be achieved in many ways,
#' mathematically most appealing by multiplication with a diagonal matrix from left
#' and right, or more efficiently by using sweep(.., FUN = "/") twice.
#' The \code{\link[stats]{cov2cor}} function is even a bit more efficient,
#' and provided mostly for didactical reasons.
#'
#' @references{
#'  \insertRef{Gower:71:00}{ProcMod}
#'
#'  \insertRef{Lingoes:74:00}{ProcMod}
#'
#'  \insertRef{Coissac-Eric:19:00}{ProcMod}
#' }
#'
#' @param ...   the set of matrices or a \code{\link[ProcMod]{procmod_frame}}
#'              object.
#' @param nrand number of randomisation used to estimate the mean
#'              covariance observed between two random matrix.
#'              If rand is \code{NULL} or equal to \code{0}, no correction
#'              is estimated and the raw procrustean covariances are
#'              estimated.
#' @param p_adjust_method the multiple test correction method used
#'              to adjust p values. \code{p_adjust_method}
#'              belongsone of the folowing values: \code{"holm"},
#'              \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'              \code{"BH"}, \code{"BY"}, \code{"fdr"},  \code{"none"}.
#'              The default is,set to \code{"holm"}.
#'
#' @return a \code{procmod_varls} object which corresponds to a numeric
#'              matrix annotated by several attributes.
#'
#'              The following attribute is always added:
#'
#'              - \code{nrand} an integer value indicating the number of
#'                randomisations used to estimate the mean of the random
#'                covariance.
#'
#'              When \code{nrand} is greater than 0 a couple of attributes
#'              is added:
#'
#'              - \code{rcovls} a numeric matrix containing the estimation
#'                of the mean of the random covariance.
#'
#'              - \code{p.value} a numeric matrix containing the estimations
#'                of the p.values of tests checking that the observed
#'                covariance is larger than the mean of the random covariance.
#'                p.values are corrected for multiple tests according to the
#'                method specified by the \code{p_adjust_method} parameter.
#'
#' @examples
#' # Build Three matrices of 3 rows.
#' A <- simulate_matrix(10,3)
#' B <- simulate_matrix(10,5)
#' C <- simulate_correlation(B,10,r2=0.6)
#'
#' # Computes the variance covariance matrix
#' varls(A, B, C)
#' varls(A = A, B = B, C = C)
#'
#' data = procmod_frame(A = A, B = B, C = C)
#' varls(data)
#'
#' # Computes the correlation matrix
#' corls(data, nrand = 500)
#'
#' # Computes the partial correlation matrix
#' corls_partial(data)
#' corls_partial(data, nrand = 0)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#'
#' @seealso \code{\link[stats]{p.adjust}}
#'
#' @rdname varls
#' @name varls
#' @aliases varls
#' @aliases corls
#' @aliases corls_partial
#' @export
varls <- function(...,
                  nrand = 100,
                  p_adjust_method = "holm") {
  if (!is.null(nrand) && length(nrand) > 1) {
    stop("nrand must be a single numeric value or the NULL value")
  }

  if (nrand == 0) nrand <- NULL

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

  nx <- length(xs)
  n <- nrow(xs)
  n1 <- n - 1

  cov_xxs <- Matrix::Matrix(0, nrow = nx, ncol = nx)

  # Computes the standard Covls covariance matrix
  for (i in seq_len(nx))
    for (j in i:nx)
      cov_xxs[i, j] <- sum(svd(crossprod(xs[[i]],xs[[j]]))$d) / n1

  cov_xxs <- as.matrix(Matrix::forceSymmetric(cov_xxs, uplo = "U"))

  # Computes Covls under null hypothesis
  if (!is.null(nrand)) {
    n_r_greater <- array(0, dim = dim(cov_xxs))

    v_xs <- vector(mode = "list", nx)
    for (i in seq_len(nx))
      v_xs[[i]] <- var(xs[[i]])

    if (.has_doParallel && getDoParRegistered()) {
        `%dp%` <- `%dopar%`
      }
    else{
        `%dp%` <- `%do%`
      }

    s_cov_xxs <- foreach(k = seq_len(nrand),
                         .combine = cbind) %dp% {
      s1_cov_xxs <- matrix(0, nrow = nx, ncol = nx)
      r_xs <- vector(mode = "list", nx)
      r_ys <- vector(mode = "list", nx)
      for (i in seq_len(nx)) {
        r_xs[[i]] <- MASS::mvrnorm(n,
                                   rep(0,nrow(v_xs[[i]])),
                                   Sigma = v_xs[[i]],
                                   empirical = TRUE)
        r_ys[[i]] <- MASS::mvrnorm(n,
                                   rep(0,nrow(v_xs[[i]])),
                                   Sigma = v_xs[[i]],
                                   empirical = TRUE)
      }

      for (i in seq_len(nx))
        for (j in i:nx)
          s1_cov_xxs[i, j] <- sum(svd(crossprod(r_xs[[i]],r_ys[[j]]))$d)

      s1_cov_xxs  / n1

                         }

    dim(s_cov_xxs) = c(nx, nx, nrand)

    s_cov_xxs <- s_cov_xxs

    r_cov_xxs <- apply(s_cov_xxs,
      MARGIN = c(1,2),
      FUN = mean
    )

    for (k in seq_len(nrand))
      n_r_greater <- n_r_greater + (s_cov_xxs[, , k] >= cov_xxs)

    r_cov_xxs <- as.matrix(Matrix::forceSymmetric(r_cov_xxs, uplo = "U"))
    cov_xxs <- pmax(cov_xxs - r_cov_xxs, 0)

    colnames(r_cov_xxs) <- x_names
    rownames(r_cov_xxs) <- x_names

    p_values <- n_r_greater / nrand

    c_p_values <- p.adjust(p_values[upper.tri(p_values, diag = FALSE)],
      method = p_adjust_method
    )

    p_values[upper.tri(p_values, diag = FALSE)] <- c_p_values
    p_values <- as.matrix(Matrix::forceSymmetric(p_values, uplo = "U"))

    colnames(p_values) <- x_names
    rownames(p_values) <- x_names

    attr(cov_xxs, "rcovls") <- as.matrix(r_cov_xxs)
    attr(cov_xxs, "p.value") <- p_values
    attr(cov_xxs, "nrand") <- nrand
  }

  colnames(cov_xxs) <- x_names
  rownames(cov_xxs) <- x_names

  make_subS3Class(cov_xxs, "procmod_varls")
}


#' @rdname varls
#' @export
corls <- function(..., nrand = 100,
                  p_adjust_method = "holm") {
  cov <- varls(...,
    nrand = nrand,
    p_adjust_method = p_adjust_method
  )
  s <- sqrt(diag(cov))
  vv <- s %o% s
  rls <- cov / vv
  class(rls) <- "matrix"


  if (!is.null(attr(cov, "rcovls"))) {
    attr(rls, "nrand") <- attr(cov, "nrand")
  }

  make_subS3Class(rls, "procmod_corls")
}

#' @rdname varls
#' @export
corls_partial <- function(..., nrand = 100) {
  rls <- corls(..., nrand = nrand)
  C <- solve(rls)
  S <- sqrt(diag(C))
  D <- S %o% S
  rp <- C / D

  attr(rp, "nrand") <- attr(rls, "nrand")
  attr(rp, "rcovls") <- attr(rls, "rcovls")

  make_subS3Class(rp, "procmod_corls")
}

#' Print procrustean Variance / Covariance Matrix.
#'
#' @param x a \code{procmod_varls}
#'          object
#' @param ... other parameters passed to other functions
#'
#' @examples
#' # Build Three matrices of 3 rows.
#' A <- simulate_matrix(10,3)
#' B <- simulate_matrix(10,5)
#' C <- simulate_correlation(B,10,r2=0.6)
#'
#' # Computes the variance covariance matrix
#' data <- procmod_frame(A = A, B = B, C = C)
#' v <- varls(data, nrand = 1000)
#'
#' print(v)
#'
#' @seealso \code{\link[ProcMod]{varls}}
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
print.procmod_varls <- function(x, ...) {
  class(x) <- "matrix"
  attr(x, "nrand") <- NULL
  attr(x, "rcovls") <- NULL
  attr(x, "p.value") <- NULL
  print(x)
}

#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
`$.procmod_varls` <- function(x, name) {
  attr(x,name)
}

#' The Names of the elements of a Variance / Covariance Matrix.
#'
#' Returns the names of the elements associated to a \code{procmod_varls}
#' object.
#'
#' @param x a \code{procmod_varls} object
#'
#' @examples
#' # Build Three matrices of 3 rows.
#' A <- simulate_matrix(10,3)
#' B <- simulate_matrix(10,5)
#' C <- simulate_correlation(B,10,r2=0.6)
#'
#' # Computes the variance covariance matrix
#' data <- procmod_frame(A = A, B = B, C = C)
#' v <- varls(data, nrand = 1000)
#'
#' names(v)
#'
#' @seealso \code{\link[ProcMod]{varls}}
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
names.procmod_varls <- function(x) {
  n <- names(attributes(x))
  bn <- grep(pattern = "^(dim|class)",
             x = n)
  n[-bn]
}


#' Print a procrustean Correlation Matrix.
#'
#' @param x a \code{procmod_corls}
#'          object
#' @param ... other parameters passed to other functions
#'
#' @examples
#' # Build Three matrices of 3 rows.
#' A <- simulate_matrix(10,3)
#' B <- simulate_matrix(10,5)
#' C <- simulate_correlation(B,10,r2=0.6)
#'
#' # Computes the correlation matrix
#' data <- procmod_frame(A = A, B = B, C = C)
#' cls <- corls(data, nrand = 1000)
#'
#' print(cls)
#'
#' @seealso \code{\link[ProcMod]{corls}}
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
print.procmod_corls <- function(x, ...) {
  class(x) <- "matrix"
  attr(x, "nrand") <- NULL
  attr(x, "rcovls") <- NULL
  attr(x, "rcorls") <- NULL
  attr(x, "p.value") <- NULL
  print(x)
}

#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
`$.procmod_corls` <- function(x, name) {
  attr(x,name)
}


#' The Names of the elements of a Correlation Matrix
#'
#' Returns the names of the elements associated to a \code{procmod_corls}
#' object.
#'
#' @param x a \code{procmod_corls} object
#'
#' @examples
#' # Build Three matrices of 3 rows.
#' A <- simulate_matrix(10,3)
#' B <- simulate_matrix(10,5)
#' C <- simulate_correlation(B,10,r2=0.6)
#'
#' # Computes the correlation matrix
#' data <- procmod_frame(A = A, B = B, C = C)
#' cls <- corls(data, nrand = 1000)
#'
#' names(cls)
#'
#' @seealso \code{\link[ProcMod]{corls}}
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
names.procmod_corls <- function(x) {
  n <- names(attributes(x))
  bn <- grep(pattern = "^(dim|class)",
             x = n)
  n[-bn]
}


