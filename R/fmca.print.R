#' @title Print \code{fmca} object
#'
#' Print a \code{fmca} object
#'
#' @param x \code{fmca} object (see \link{compute_optimal_encoding} function)
#' @param n maximal number of rows and cols to print
#' @param ... Not used.
#'
#'
#' @method print fmca
#'
#' @seealso \link{compute_optimal_encoding} \link{summary.fmca}
#'
#' @export
print.fmca <- function(x, n = 6, ...) {
  cat("$basisobj:\n")
  cat("  Type:", x$basisobj$type, "\n")
  cat("  Range:", x$basisobj$rangeval[1], "to", x$basisobj$rangeval[2], "\n")
  cat("  Number of basis functions:", x$basisobj$nbasis, "\n")
  cat("\n$label:", nrow(x$label), "elements", "\n")
  print(head(x$label$label, n))
  cat("\n$eigenvalues:", length(x$eigenvalues), "elements", "\n")
  print(head(x$eigenvalues, n))
  cat("\n$alpha:", "list of", length(x$alpha), "elements", "\n")
  cat("$alpha[[1]]:", nrow(x$alpha[[1]]), "rows", ncol(x$alpha[[1]]), "columns", "\n")
  print(head(x$alpha[[1]], n))
  cat("\n$pc:", nrow(x$pc), "rows", ncol(x$pc), "columns", "\n")
  print(head(x$pc[, 1:n], n))
  cat("\nOther elements: \"F\", \"G\", \"V\", \"pt\"\n")
}


#' @title Object Summaries
#'
#' Summary of a \code{fmca} object
#'
#'
#' @param object \code{fmca} object (see \link{compute_optimal_encoding} function)
#' @param n maximal number of rows and cols to print
#' @param ... Not used.
#'
#'
#' @method summary fmca
#'
#' @seealso \link{compute_optimal_encoding} \link{print.fmca}
#'
#' @export
summary.fmca <- function(object, n = 6, ...) {
  cat("#### FMCA\n\n")
  cat("## Data \n")
  cat("Number of individuals:", nrow(object$pc), "\n")
  cat("Number of states:", nrow(object$label), "\n")
  cat("Time Range:", object$basisobj$rangeval[1], "to", object$basisobj$rangeval[2], "\n")
  cat("States:  ")
  cat(head(object$label$label, n), "\n")
  cat("\n")
  cat("## Basis \n")
  cat("Type:", object$basisobj$type, "\n")
  cat("Number of basis functions:", object$basisobj$nbasis, "\n")
  cat("\n")
  cat("## Outputs\n")
  cat("Eigenvalues:\n  ")
  cat(head(object$eigenvalues, n), "\n")
  cat("\nExplained variance:\n  ")
  cat(round(head(cumsum(object$eigenvalues)/sum(object$eigenvalues), n), 3), "\n")
  cat("\nOptimal encoding:\n")
  print(head(object$alpha[[1]], n))
  cat("\nPrincipal components:\n")
  print(head(object$pc[, 1:n]))
  cat("\n")
  cat("Total elapsed time:", round(object$runTime, 3), "s\n")
}
