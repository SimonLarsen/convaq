#' Perform CNV-based association study.
#' @export
#' @param segments1 Data frame of segments for group 1.
#' @param segments2 Data frame of segments for group 2.
#' @param model Model type. Either "statistical" or "query".
#' @param qvalues TRUE if q-values should be computed, FALSE otherwise.
#' @param qvalues.rep Number of repetitions to use in q-value computation.
#' @param nthreads Number of threads to use. Defaults to number of cores available.
#' @param p.cutoff (statistical) P-value cutoff in statistical model.
#' @param pred1 (query) Predicate for group 1 in query model.
#' @param pred2 (query) Predicate for group 2 in query model.
convaq <- function(
  segments1,
  segments2,
  model,
  name1 = "Group 1",
  name2 = "Group 2",
  qvalues = FALSE,
  qvalues.rep = 500,
  nthreads = NULL,
  p.cutoff = 0.05,
  pred1 = NULL,
  pred2 = NULL
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
  
  # convert model to numeric value
  model.num <- match(model.full, c("statistical","query"))

  if(is.null(nthreads)) nthreads <- 0
  
  comp1 <- 0; comp2 <- 0; value1 <- 0; value2 <- 0
  eq1 <- 0; eq2 <- 0; type1 <- 0; type2 <- 0

  if(model.full == "query") {
    if(is.null(pred1)) stop("Missing predicate for group 1.")
    if(is.null(pred2)) stop("Missing predicate for group 2.")
    pred1_l <- parse_predicate(pred1, types)
    pred2_l <- parse_predicate(pred2, types)
    comp1 <- pred1_l$comp; value1 <- pred1_l$value; eq1 <- pred1_l$eq; type1 <- pred1_l$type
    comp2 <- pred2_l$comp; value2 <- pred2_l$value; eq2 <- pred2_l$eq; type2 <- pred2_l$type
  }
  
  # call C++ backend
  out <- convaqCpp(
    segments1, segments2,
    model.num,
    qvalues, qvalues.rep,
    nthreads,
    p.cutoff,
    comp1, value1, eq1, type1,
    comp2, value2, eq2, type2
  );
  
  if(is.null(out) || length(out) == 0) {
    message("No variations found.")
    return(NA)
  }
  
  out <- data.frame(out)
  colnames(out) <-  c(
    "chr", "start", "end", "length", "type", "pvalue", "qvalue",
    paste0(name1, ": ", types.pretty),
    paste0(name2, ": ", types.pretty)
  )

  # convert
  out$type <- factor(types.pretty[out$type+1], levels=c(types.pretty,"Normal"))

  # set p-values to NA if using query model
  if(model.full == "query") {
    out$pvalue <- NA
  }
  # set q-values to NA if not computed
  if(!qvalues) {
    out$qvalue <- NA
  }

  return(out)
}
