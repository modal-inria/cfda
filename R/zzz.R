#' @export
.onAttach <- function(lib, pkg) {
  packageStartupMessage(
    paste(
      "INFO: Since cfda 0.11, the computing method of compute_optimal_encoding function has changed.",
      "To use the former one, add the argument: method='parallel'"
    )
  )
}