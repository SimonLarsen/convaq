#' @export
states <- function(x, ...) UseMethod("states")

#' Extract states of individual samples from CoNVaQ result as a table.
#' 
#' @param x A convaq object.
#' @export
states.convaq <- function(x) {
  if(class(x) != "convaq") stop("Object is not a convaq object.")
  
  fix <- function(x) {
    if(length(x) == 0) NA else paste0(x, collapse=",")
  }
  
  rows <- lapply(x$state, function(re) {
    r1 <- rbind(sapply(re[[1]], fix))
    r2 <- rbind(sapply(re[[2]], fix))
    cbind(r1, r2)
  })
  
  out <- do.call(rbind, rows)
  data.frame(out)
}
