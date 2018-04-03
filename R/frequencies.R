#' Extract variation frequencies for each region as a table.
#' 
#' @param x A convaq object.
#' @param ... Further arguments passed to or from other methods.
#' @export
frequencies <- function(x, ...) UseMethod("frequencies")

#' Extract variation frequencies for each region as a table.
#' 
#' @param x A convaq object.
#' @param ... Further arguments passed to or from other methods.
#' @export
frequencies.convaq <- function(x, ...) {
  if(class(x) != "convaq") stop("Object is not a convaq object.")
  
  fix <- function(s) {
    su <- unique(signif(unlist(s)*100, digits=3))
    if(length(su) == 1) sprintf("%s %%", formatC(su))
    else sprintf("%s-%s %%", formatC(min(su)), formatC(max(su)))
  }
  
  types.pretty <- c("Gain","Loss","LOH")
  rows <- lapply(x$freq, function(re) {
    r1 <- rbind(sapply(re[[1]], fix))
    r2 <- rbind(sapply(re[[2]], fix))
    colnames(r1) <- paste0(x$name1, ": ", types.pretty)
    colnames(r2) <- paste0(x$name2, ": ", types.pretty)
    data.frame(r1, r2, check.names=FALSE)
  })
  
  data.frame(do.call(rbind, rows), check.names=FALSE)
}
