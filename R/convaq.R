#' Perform CNV-based association study.
#' @export
convaq <- function(
  segments1,
  segments2,
  model,
  p.cutoff = 0.05,
  qvalues = FALSE,
  qvalues.rep = 500,
  nthreads = NULL
) {
  # convert segment types to numbers.
  # Gain = 0, Loss = 1, LOH = 2.
  types.pretty <- c("Gain","Loss","LOH")
  types <- tolower(types.pretty)

  segments1$type <- tolower(segments1$type)
  segments2$type <- tolower(segments2$type)

  found.types <- unique(c(segments1$type, segments2$type))
  bad.types <- found.types[!(found.types %in% types)]
  if(length(bad.types) > 0) {
    stop("Invalid segment type(s): ", paste0(bad.types, collapse=", "))
  }

  segments1$type <- as.numeric(as.factor(segments1$type))-1
  segments2$type <- as.numeric(as.factor(segments2$type))-1

  # convert patients to numbers 0, 1, ...
  patients1 <- as.factor(segments1$patient)
  segments1$patient <- as.numeric(patients1)-1
  patients2 <- as.factor(segments2$patient)
  segments2$patient <- as.numeric(patients2)-1

  # convert chromosomes to strings
  segments1$chr <- as.character(segments1$chr)
  segments2$chr <- as.character(segments2$chr)

  model.full <- tryCatch(
    match.arg(model, c("statistical","query")),
    error = function(e) NULL
  )
  if(is.null(model.full)) stop("Unrecognized model type: ", model)

  # convert model to number
  model.num <- switch(model.full,
    statistical = 0,
    query = 1
  )

  if(is.null(nthreads)) nthreads <- 0

  # call C++ backend
  out <- convaqCpp(segments1, segments2, model.num, p.cutoff, qvalues, qvalues.rep, nthreads);

  # convert
  out$type <- factor(types.pretty[out$type+1], levels=types.pretty)

  # set q-values to NA if not computed
  if(!qvalues) {
    out$qvalue <- NA
  }

  return(out)
}
