#' @include internals.R
NULL

#' Internal function repeating a matrix.
#'
#' @description repeats several times the rows of a matrix
#'              to create a new matrix with more rows. The
#'              final row count must be a multiple of the
#'              initial row count
#'
#' @param x  The matrix to replicate
#' @param nrow  an interger value specifying the number of row
#'                    of the returned matrix
#'
#' @return a new matrix with the same number of columns but with `nrow`
#'         rows.
#' @note Internal function do not use.
#'
#' @author Eric Coissac <eric.coissac@metabarcoding.org>
#' @author Christelle Gonindard-Melodelima <christelle.gonindard@metabarcoding.org>
#' @rdname internal.rep_matrix
#'
.rep_matrix <- function(x, nrow) {
  N <- nrow(x)

  if ((nrow %% N != 0L)) {
    stop(sprintf(
      "The size of the longest object (%d) is not a multiple of the size of the shortest (%d)",
      nrow, N
    ),
    domain = NA
    )
  }

  rep <- x
  while (nrow(rep) < nrow)
    rep <- rbind(rep, x)

  return(rep)
}

#' Internal function coercing the data to a matrix.
#'
#' @description Transforme the \code{x} value into a \code{numeric matrix} of
#'              the correct size or into a \code{dist} object.
#'
#' @param x  The data to coerce
#' @param nrows  an interger value specifying the number of row
#'               of the returned matrix
#' @param contrasts see the \code{contrasts_arg} argument
#'               of the \code{\link[ProcMod]{procmod_frame}}
#'               constructor.
#'
#' @return a new numeric matrix with correct size.
#' @note Internal function do not use.
#'
#' @author Eric Coissac <eric.coissac@metabarcoding.org>
#' @author Christelle Gonindard-Melodelima <christelle.gonindard@metabarcoding.org>
#' @rdname internal.procmod_coerce_value
#'
.procmod_coerce_value <- function(x, nrows = 0, contrasts = NULL) {
  xi <- if (is.data.frame(x)) {
    as.matrix(x)
  } else if (is.matrix(x) || inherits(x, "dist")) {
    x
  } else if (is.factor(x)) {
    if (is.null(contrasts)) {
      contrasts(x)[x, ]
    } else if (is.character(contrasts) || is.function(contrasts)) {
      match.fun(contrasts, descend = FALSE)(x)[x, ]
    } else {
      contrasts[x, ]
    }
  }
  else {
    as.matrix(x)
  }

  if (!(is.matrix(xi) || inherits(x, "dist"))) {
    stop("Value cannot be coerced to a Matrix")
  }

  if (is.matrix(xi)) {
    dxi <- dim(xi)
    rownamesi <- rownames(xi)
    colnamesi <- colnames(xi)

    xi <- as.numeric(xi)

    dim(xi) <- dxi
    rownames(xi) <- rownamesi
    colnames(xi) <- colnamesi

    N <- nrow(xi)

    if (nrows > 0L && N < nrows) {
      if (N > 0L && (nrows %% N == 0L)) {
        xi <- .rep_matrix(xi, nrow = nrows)
      } else {
        stop(sprintf(
          ngettext(
            N, "replacement has %d row, data has %d",
            "replacement has %d rows, data has %d"
          ),
          N, nrows
        ),
        domain = NA
        )
      }
    }
  }
  else {
    N <- attr(xi, "Size")
    if (nrows > 0L && N != nrows) {
      stop(sprintf(
        ngettext(
          N, "replacement has %d row, data has %d",
          "replacement has %d rows, data has %d"
        ),
        N, nrows
      ),
      domain = NA
      )
    }
  }

  return(xi)
}

.siteNames <- function(x) {
  if (inherits(x, "dist")) {
    attr(x, "Labels")
  } else {
    rownames(x)
  }
}

`.siteNames<-` <- function(x, value) {
  if (inherits(x, "dist")) {
    stopifnot(is.null(value) || length(value) == attr(x, "Size"))
    attr(x, "Labels") <- value
  }
  else {
    rownames(x) <- value
  }

  x
}


.siteCount <- function(x) {
  if (inherits(x, "dist")) {
    attr(x, "Size")
  } else {
    nrow(x)
  }
}


.siteSelect <- function(x, select) {
  if (inherits(x, "dist")) {
    as.dist(as.matrix(x)[select, select, drop = FALSE])
  } else {
    x[select, , drop = FALSE]
  }
}


#' The procmod_frame data structure.
#'
#' A \code{procmod_frame} can be considered as the analog of a
#' \code{data.frame} for vector data. In a \code{procmod_frame}
#' each element, equivalent to a column in a \code{data.frame}
#' is a numeric matrix or a distance matrix object (\code{dist}).
#' Every element must describe the same number of individuals.
#' Therefore every numeric matrix must have the same number of row
#' (\code{nrow}) and every distance matrix must have the same size
#' (\code{attr(d,"Size")}). A \code{procmod_frame} can simultaneously
#' contain both types of data, numeric and distance matrix.
#'
#' @param ... a set of objects to aggregate into a
#'        \code{procmod_frame}. These objects can be
#'        numeric matrices, or dist objects. Every objects
#'        must have the same number of row.
#'
#' @param row_names a character vector containing names associated
#'        to each row.
#'
#' @param check_rows a logical value. When set to \code{TRUE}, its
#'        default value, the number of row of every elements of the
#'        \code{procmod_frame} are tested for equality. Otherwise no
#'        check is done.
#'
#' @param reorder_rows a logical value. When set to \code{TRUE}, its
#'        default value, every elements of the
#'        \code{procmod_frame} are reordered according to the \code{row_names}
#'        order. Otherwise nothing is done.
#'
#' @param contrasts_arg	A list, whose entries are values
#'        (numeric matrices or character strings naming functions)
#'        to be used as replacement values for the contrasts
#'        replacement function and whose names are the names
#'        of columns of data containing factors.
#'
#' @return a \code{procmod_frame} instance.
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
#' length(dataset)
#' nrow(dataset)
#' ncol(dataset)
#' dataset$euk
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
procmod_frame <- function(...,
                          row_names = NULL,
                          check_rows = TRUE,
                          reorder_rows = TRUE,
                          contrasts_arg = NULL) {
  has.row_names <- !missing(row_names)

  varnames <- dots_names(...)

  x <- list(...)
  n <- length(x)

  #  message(x)
  #  message(n)
  #  message(row_names)

  if ((!has.row_names || is.null(row_names)) && n >= 1) {
    row_names <- .siteNames(x[[1]])
  }

  nrows <- integer(n)
  value <- vector(mode = "list", length = n)

  names(value) <- varnames

  # types <- character(n)

  for (i in seq_len(n)) {
    contrasts <- contrasts_arg[varnames[i]]
    if (i == 1) {
      xi <- .procmod_coerce_value(x[[i]],
        contrasts = contrasts
      )
    } else {
      xi <- .procmod_coerce_value(x[[i]],
        nrows = nrows[1],
        contrasts = contrasts
      )
    }

    if (reorder_rows &&
      !is.null(row_names) &&
      !is.null(.siteNames(xi))) {
      xi <- .siteSelect(xi, row_names)
    }

    nrows[i] <- .siteCount(xi)
    value[[i]] <- xi
  }

  stopifnot(all(nrows[i] == nrows))
  # message(row_names," : ",length(row_names),",",nrows[i])
  if (length(row_names) == nrows[i]) {
    attr(value, "row_names") <- row_names
    if (check_rows) {
      for (i in seq_len(n)) {
        if (!all(row_names == .siteNames(value[[i]]))) {
          stop("Row names among matrices are not consistant")
        }
      }
    } else {
      for (i in seq_len(n))
        .siteNames(value[[i]]) <- row_names
    }
  }
  else {
    for (i in seq_len(n))
      .siteNames(value[[i]]) <- NULL
  }

  make_subS3Class(value, "procmod_frame")
}

#'
#' Check if an object is a ProcMod Frame.
#'
#' @param x a R \code{object to test}
#'
#' @return a \code{logical} value equals to
#'         \code{TRUE} if \code{x} is a \code{procmod_frame},
#'         \code{FALSE} otherwise.
#'
#' @examples
#' # Builds a procmod_frame with two random matrices
#' m1 <- simulate_matrix(10,20)
#' m2 <- simulate_matrix(10,30)
#' pmf <- procmod_frame(m1 = m1, m2 = m2)
#'
#' # Returns TRUE
#' is_procmod_frame(pmf)
#'
#' # Returns FALSE
#' is_procmod_frame(3)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
is_procmod_frame <- function(x) {
  inherits(x, "procmod_frame")
}

#'
#' Coerce to a ProcMod Frame.
#'
#' Conversion methods are proposed for \code{list},
#' \code{matrix} and \code{array}.
#'
#' @param data a R object to coerce.
#' @param ... supplementary parameters used in some
#'            implementation of that method
#' @return a \code{procmod_frame} object
#'
#' @examples
#' # Builds a list containing two random matrices
#' m1 <- simulate_matrix(10,20)
#' m2 <- simulate_matrix(10,30)
#' l <- list(m1 = m1, m2 = m2)
#'
#' # Converts the list to a procmod_frame
#' pmf1 <- as_procmod_frame(l)
#'
#' # Builds a procmod_frame from a matrix
#' m3 <- matrix(1:12,nrow=3)
#' pmf2 <- as_procmod_frame(matrix(1:12,nrow=3))
#' # Returns 4, the column count of the input matrix
#' length(pmf2)
#'
#' # Builds a 3D array
#' a <- array(1:24,dim = c(3,4,2))
#'
#' # The conversion to a procmod_frame makes
#' # an procmod element from each third dimension
#' as_procmod_frame(a)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
as_procmod_frame <- function(data, ...) {
  UseMethod("as_procmod_frame", data)
}

#' @rdname as_procmod_frame
#' @export
as_procmod_frame.list <- function(data, ...) {
  data <- c(data, list(...))
  do.call(procmod_frame, data)
}

#' @rdname as_procmod_frame
#' @export
as_procmod_frame.procmod_frame <- function(data, ...) {
  data
}

#' @rdname as_procmod_frame
#' @export
as_procmod_frame.array <- function(data, ...) {
  di <- dim(data)
  stopifnot(length(di) == 3)

  l <- lapply(
    seq_len(di[3]),
    function(i) data[, , i]
  )

  if (length(attr(data, "dimnames")) == 3) {
    names(l) <- attr(data, "dimnames")[[3]]
  }

  data <- c(l, list(...))

  do.call(procmod_frame, l)
}

# #' @rdname procmod_frame
# #' @export
# as_procmod_frame.pm <- function(data, ...) {
#   vars.procmod(terms(data), data$model)
# }

#' @rdname as_procmod_frame
#' @export
as_procmod_frame.matrix <- function(data, ...) {
  l <- vector(mode = "list", length = ncol(data))
  for (i in seq_len(ncol(data))) {
    l[[i]] <- data[, i]
  }

  if (!is.null(colnames(data))) {
    names(l) <- colnames(data)
  }

  as_procmod_frame(l)
}

#' Dimensions of a ProcMod Frame.
#'
#' Dimension 1 is the number of rows (individus)
#' shared by the aggregated matrices. Dimension 2
#' is the number of aggregated matrices
#'
#' @param x a \code{\link[ProcMod]{procmod_frame}}
#'          object
#'
#' @examples
#' # Builds a procmod_frame with two random matrices
#' m1 <- simulate_matrix(10,20)
#' m2 <- simulate_matrix(10,30)
#' pmf <- procmod_frame(m1 = m1, m2 = m2)
#' dim(pmf)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
dim.procmod_frame <- function(x)
  return(c(.siteCount(x[[1]]), length(x)))

#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
`[[<-.procmod_frame` <- function(x, i, value) {
  cl <- class(x)
  nrows <- .siteCount(x)
  class(x) <- "list"

  if (!is.null(value)) {
    value <- .procmod_coerce_value(value, nrows)
    N <- .siteCount(value)

    if (N > nrows) {
      stop(sprintf(
        ngettext(
          N, "replacement has %d row, data has %d",
          "replacement has %d rows, data has %d"
        ),
        N, nrows
      ),
      domain = NA
      )
    }

    if (N < nrows) {
      stop(sprintf(
        ngettext(
          N, "replacement has %d row, data has %d",
          "replacement has %d rows, data has %d"
        ),
        N, nrows
      ),
      domain = NA
      )
    }
  }

  .siteNames(value) <- attr(x, "row_names")

  x[[i]] <- value
  class(x) <- cl

  return(x)
}

#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
`$<-.procmod_frame` <- function(x, name, value) {
  x[[name]] <- value

  x
}

#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
`[.procmod_frame` <- function(x, i, j,
                              drop = TRUE) {
  has.j <- !missing(j)
  has.i <- !missing(i)
  has.drop <- !missing(drop)
  Narg <- nargs() - 2 + (has.i | has.j | has.drop) - has.drop

  # message("Nargs = ",Narg," i:",has.i," j:",has.j," drop:",has.drop)

  if (!all(names(sys.call()) %in% c("", "drop"))) {
    warning("named arguments other than 'drop' are discouraged")
  }


  if (!has.i && !has.j) {
    # Case 1 : X[]
    return(x)
  } else if (!has.i && has.j) {
    # Case 2 : X[,j] ou x[i]
    # message('Case 2 : X[,j]')
    y <- as.list(x)[j]
  }
  else if (has.i && Narg == 1) {
    # Case 3 : X[,j] ou x[i]
    # message('Case 3 : X[i]')
    y <- as.list(x)[i]
  }
  else if (has.i && !has.j && Narg > 1) {
    # Case 4 : X[i,]
    # message('Case 4 : X[i,]')
    y <- lapply(x, function(m) .siteSelect(m, i))
  }
  else if (has.i && has.j) {
    # message('Case 5 : X[i,j]')
    y <- x[j, drop = FALSE]
    y <- y[i, , drop = FALSE]
  }

  if (drop && length(y) == 1L) {
    y <- y[[1]]
  } else {
    y <- make_subS3Class(y, "procmod_frame")
    attr(y, "row_names") <- rownames(y[[1]])
  }

  y
}

#' Subsetting Procmod Frames
#'
#' This is the implementation of the \code{\link[base]{subset}} generic function for
#' \code{procmod_frame}.
#'
#' The subset argument works on rows. Note that subset will be evaluated in the
#' \code{procmod_frame}, so columns can be referred to (by name) as variables
#' in the expression (see the examples).
#'
#' The select argument if provided indicates with matrices
#' have to be conserved.  It works by first replacing column names in the selection
#' expression with the corresponding column numbers in the \code{procmod_frame} and then using
#' the resulting integer vector to index the columns. This allows the use of the
#' standard indexing conventions so that for example ranges of columns can
#' be specified easily, or single columns can be dropped (see the examples). Remember
#' that each column of a \code{procmod_frame} is actually a matrix.
#'
#' The drop argument is passed on to the \code{procmod_frame} indexing method.
#' The default value is \code{FALSE}.
#'
#' @param x	object to be subsetted.
#' @param subset	logical expression indicating elements or
#'                rows to keep: missing values are taken as false.
#' @param select	expression, indicating columns to select from a data frame.
#' @param drop	  passed on to [ indexing operator.
#' @param ...	further arguments to be passed to or from other methods.
#'
#' @return A \code{procmod_frame} containing just the selected rows and columns.
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
#' dim(dataset)
#'
#' higher_ph = subset(dataset,soil[,"pH"] > 0)
#' dim(higher_ph)
#'
#' without_bacteria = subset(dataset,soil[,"pH"] > 0, -bac)
#' dim(without_bacteria)
#'
#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
subset.procmod_frame <- function(x, subset, select, drop = FALSE, ...) {
  r <- if (missing(subset)) {
    rep_len(TRUE, nrow(x))
  } else {
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r)) {
      stop("'subset' must be logical")
    }
    r & !is.na(r)
  }
  vars <- if (missing(select)) {
    TRUE
  } else {
    nl <- as.list(seq_along(x))
    names(nl) <- names(x)
    eval(substitute(select), nl, parent.frame())
  }
  x[r, vars, drop = drop]
}

#' @author Eric Coissac
#' @author Christelle Gonindard-Melodelima
#' @export
as.list.procmod_frame <- function(x, ...) {
  class(x) <- "list"

  x
}
