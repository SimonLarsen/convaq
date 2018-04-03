#' Extract states of individual samples from CoNVaQ result as a table.
#' 
#' @param x A convaq object.
#' @param ... Further arguments passed to or from other methods.
#' @export
states <- function(x, ...) UseMethod("states")

#' Extract states of individual samples from CoNVaQ result as a table.
#' 
#' @param x A convaq object.
#' @param ... Further arguments passed to or from other methods.
#' @export
states.convaq <- function(x, ...) {
  if(class(x) != "convaq") stop("Object is not a convaq object.")
  
  rows <- lapply(x$state, function(re) {
    r1 <- rbind(sapply(re[[1]], paste0, collapse=","))
    r2 <- rbind(sapply(re[[2]], paste0, collapse=","))
    cbind(r1, r2)
  })
  
  data.frame(do.call(rbind, rows))
}
