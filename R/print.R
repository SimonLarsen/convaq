#' Print description of CoNVaQ results object.
#'
#' @param x A convaq object.
#' @export
print.convaq <- function(x) {
  cat("CoNVaQ results object.\n\n")
  cat("Model:               ", x$model, "\n")
  cat("Group 1 name:        ", x$name1, "\n")
  cat("Group 2 name:        ", x$name2, "\n")
  cat("No. regions found:   ", nrow(x$regions), "\n")
  cat("Compute q-values:    ", x$qvalues, "\n")
  cat("Q-value repetitions: ", x$qvalues.rep, "\n")
  if(x$model == "statistical") {
  cat("P-value cutoff:      ", x$p.cutoff, "\n")
  }
  else if(x$model == "query") {
  cat("Predicate 1:         ", x$pred1, "\n")
  cat("Predicate 2:         ", x$pred2, "\n")
  }
}
