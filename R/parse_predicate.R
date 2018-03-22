parse_predicate <- function(pred, types) {
  parts <- unlist(strsplit(tolower(pred), " "))
  if(length(parts) != 4) stop("Malformed predicate.")
  
  comp <- match(parts[1], c("<",">","<=",">="))
  value <- as.numeric(parts[2])
  eq <- match(parts[3], c("==","!="))
  type <- match(parts[4], c(types, "normal"))-1
  
  if(is.na(comp)) stop("Invalid operator in predicate: ", parts[1])
  if(is.na(value)) stop("Invalid value in predicate: ", parts[2])
  if(value < 0 | value > 1) stop("Invalid value in predicate: ", value, ". Must be between 0 and 1.")
  if(is.na(type)) stop("Invalid variation type in predicate: ", type)
  
  list(comp=comp, value=value, eq=eq, type=type)
}