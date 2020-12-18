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
#' @seealso \link{compute_optimal_encoding}
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
  cat("\nOther elements: \"F\", \"G\", \"V\", \"pt\"")
}
