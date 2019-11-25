make_subS3Class <- function(obj, subclass) {
  class(obj) <- c(
    paste(subclass,
      collapse = "_"
    ),
    class(obj)
  )
  return(obj)
}

dots_names <- function(...) {
  varnames <- substitute(list(...))[-1L]
  dots <- list(...)
  isname <- sapply(varnames, is.name)
  charname <- as.character(varnames)
  charname[!isname] <- ""

  n <- length(dots)

  explicit <- names(dots)

  if (is.null(explicit)) {
    explicit <- character(n)
  }

  ze <- !nzchar(explicit)

  explicit[ze] <- charname[ze]
  ze <- !nzchar(explicit)

  dnames <- paste("V", seq_len(n), sep = "")
  explicit[ze] <- dnames[ze]

  return(explicit)
}

make_procmod_subS3Class <- function(obj, subclass) {
  class(obj) <- c(
    paste("procmod", subclass,
      sep = "_", collapse = "_"
    ),
    class(obj)
  )

  return(obj)
}

make_procmod_data <- function(obj, subclass) {
  eud <- inherits(obj, "procmod_data", which = TRUE)

  if (eud > 0) {
    class(obj) <- class(obj)[-1:-(eud - 1)]
  } else {
    obj <- make_procmod_subS3Class(obj, "data")
  }

  if (!missing(subclass)) {
    obj <- make_procmod_subS3Class(obj, subclass)
  }

  return(obj)
}
