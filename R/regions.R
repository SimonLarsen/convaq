#' @export
regions <- function(x, ...) UseMethod("regions")

#' Extract matching CNV regions from CoNVaQ result as a table.
#'
#' @param x A convaq object.
#' @export
regions.convaq <- function(x) {
  if(class(x) != "convaq") stop("Object is not a convaq object.")

  x$regions
}
