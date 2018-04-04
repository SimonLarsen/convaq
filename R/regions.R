#' Extract matching CNV regions from CoNVaQ result as a data frame.
#'
#' @param x A convaq object.
#' @param ... Further arguments passed to or from other methods.
#' @export
regions <- function(x, ...) UseMethod("regions")

#' Extract matching CNV regions from CoNVaQ result as a data frame.
#'
#' @param x A convaq object.
#' @param ... Further arguments passed to or from other methods.
#' @export
regions.convaq <- function(x, ...) {
  if(class(x) != "convaq") stop("Object is not a convaq object.")

  x$regions
}
