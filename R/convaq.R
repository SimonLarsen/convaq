#' Perform CNV-based association study.
#' 
#' CoNVaQ is a method for performing CNV-based association studies. It provides two models:
#' a query-based model and a standard statistical model using Fisher's exact test.
#' 
#' @section Segment data format:
#' The CNV segment sets \code{segments1} and \code{segments2} must be data frame objects with the following five columns:
#' \describe{
#'   \item{patient}{Patient identifier for the patient/sample the segment was found in.}
#'   \item{chr}{Chromosome the segment is located in.}
#'   \item{start}{First position of the segment in base pairs.}
#'   \item{end}{Last position of the segment in base pairs.}
#'   \item{type}{Segment type. One of "Gain", "Loss" or "LOH".}
#' }
#' The order of the columns is used in order to determine the contents, not the column names.
#' 
#' @section Predicate format:
#' Predicates are given as a character vector with the following format:
#' \preformatted{"[COMP] [FREQ] [EQ] [TYPE]"}
#' where:
#' \describe{
#'   \item{COMP}{is a comparison operator. One of "<" (less than), ">" (greater than), "<=" (less than or equal to) or ">=" (greater than or equal to).}
#'   \item{FREQ}{is a numerical value between 0 and 1.}
#'   \item{EQ}{is either "==" (equal to) or "!=" (not equal to).}
#'   \item{TYPE}{is a segment type. One of "Gain", "Loss", "LOH" or "Normal".}
#' }
#' The four arguments must be separated by a single space.
#' 
#' A predicate is evaluated for a specific group in a each region, and will either evaluate to TRUE or FALSE.
#' If we want to find regions where at least 50\% of the patients in a group are reported as "Gain", we can use the predicate:
#' \preformatted{">= 0.5 == Gain"}
#' Likewise if we are searching for regions where less than 25\% of the patients in a group have any kind of variation,
#' we can use the predicate
#' \preformatted{"< 0.25 != Normal"}
#' The query model works by combining two predicates, one for each group of segments.
#' For instance, we can combine two predicates to search for regions where at least 50\% of patients in the first group have a "Gain",
#' and less than 25\% of patients in the second group:
#' \preformatted{convaq(s1, s2, model="query", pred1=">= 0.5 == Gain", pred2="< 0.25 == Gain")}
#' 
#' @examples
#' data("example", package="convaq")
#' s1 <- example$disease
#' s2 <- example$healthy
#' 
#' # statistical model
#' convaq(s1, s2, model="statistical", p.cutoff=0.05, qvalues=FALSE)
#' convaq(s1, s2, model="statistical", p.cutoff=0.05, qvalues=TRUE, qvalues.rep=2000)
#' 
#' # query model
#' convaq(s1, s2, model="query", pred1=">= 0.5 == Gain", pred2="<= 0.2 == Gain")
#' convaq(s1, s2, model="query", pred1=">= 0.6 != Normal", pred2=">= 0.6 == Normal")
#' 
#' @param segments1 Data frame of segments for group 1. See details.
#' @param segments2 Data frame of segments for group 2. See details.
#' @param model Model type. Either "statistical" or "query".
#' @param name1 Name of first group.
#' @param name2 Name of second group.
#' @param qvalues TRUE if q-values should be computed, FALSE otherwise.
#' @param qvalues.rep Number of repetitions to use in q-value computation.
#' @param nthreads Number of threads to use. Defaults to number of cores available.
#' @param p.cutoff (statistical model) P-value cutoff in statistical model.
#' @param pred1 (query model) Predicate for group 1 in query model.
#' @param pred2 (query model) Predicate for group 2 in query model.
#' @return An object of class \code{convaq} with the following elements:
#'   \item{regions}{Data frame of significant regions.}
#'   \item{freq}{Data frame of within-group variation frequencies for each reported region.}
#'   \item{state}{The states of individual patients/samples for each region.}
#'   \item{model}{Type of model used.}
#'   \item{name1}{Name of first group.}
#'   \item{name2}{Name of second group.}
#'   \item{qvalues}{True if q-values were computed.}
#'   \item{qvalues.rep}{Number of repetitions used in q-value computation.}
#'   \item{merge}{True if adjacent regions of same type should be merged.}
#'   \item{merge.threshold}{Maximum distance (in base pairs) allowed between merged regions.}
#'   \item{p.cutoff}{P-value cutoff (statistical model only).}
#'   \item{pred1}{Predicate for group 1 (query model only).}
#'   \item{pred2}{Predicate for group 2 (query model only).}
#' @export
convaq <- function(
  segments1,
  segments2,
  model,
  name1 = "Group 1",
  name2 = "Group 2",
  qvalues = FALSE,
  qvalues.rep = 500,
  nthreads = NULL,
  merge = FALSE,
  merge.threshold = 0,
  p.cutoff = 0.05,
  pred1 = NULL,
  pred2 = NULL
) {
  # convert segment types to numbers.
  # Gain = 0, Loss = 1, LOH = 2.
  types.pretty <- c("Gain","Loss","LOH")
  types <- tolower(types.pretty)
  # add 3 = "Normal".
  types.pretty.full <- c(types.pretty, "Normal")
  
  # check valid number of columns
  if(ncol(segments1) < 5) stop("segments1 does not have 5 columns")
  if(ncol(segments2) < 5) stop("segments2 does not have 5 columns")
  
  # check group names are not the same
  if(name1 == name2) {
    stop("Group names cannot be identifical.")
  }
  
  # extract first five columns and set colnames
  segments1 <- segments1[,1:5]
  segments2 <- segments2[,1:5]
  colnames(segments1) <- c("patient","chr","start","end","type")
  colnames(segments2) <- c("patient","chr","start","end","type")
  
  # sanitize segment types and check validity
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
  
  comp1 <- 0; value1 <- 0; eq1 <- 0; type1 <- 0;
  comp2 <- 0; value2 <- 0; eq2 <- 0; type2 <- 0;

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
    merge, merge.threshold,
    nthreads,
    p.cutoff,
    comp1, value1, eq1, type1,
    comp2, value2, eq2, type2
  );
  
  if(is.null(out)) {
    message("No variations found.")
    return(NA)
  }

  # convert
  out$regions$type <- factor(types.pretty[out$regions$type+1], levels=c(types.pretty,"Normal"))

  #remove unnecessary columns
  remove.cols <- c()
  if(model.full == "query") {
    remove.cols <- c(remove.cols, c("pvalue", "type"))
  }
  if(!qvalues) {
    remove.cols <- c(remove.cols, "qvalue")
  }
  if(length(remove.cols) > 0) {
    out$regions <- out$regions[,-match(remove.cols, colnames(out$regions))]
  }
  
  # set names for freq object
  for(i in 1:length(out$freq)) {
    names(out$freq[[i]]) <- c(name1, name2)
    names(out$freq[[i]][[1]]) <- types.pretty
    names(out$freq[[i]][[2]]) <- types.pretty
  }
  
  # set names for state object
  for(i in 1:length(out$freq)) {
    names(out$state[[i]]) <- c(name1, name2)
    names(out$state[[i]][[1]]) <- levels(patients1)
    names(out$state[[i]][[2]]) <- levels(patients2)
    out$state[[i]][[1]] <- lapply(out$state[[i]][[1]], function(x) types.pretty.full[x+1])
    out$state[[i]][[2]] <- lapply(out$state[[i]][[2]], function(x) types.pretty.full[x+1])
  }
  
  # create output object
  result <- list()
  result$model <- model.full
  result$name1 <- name1
  result$name2 <- name2
  result$regions <- out$regions
  result$freq <- out$freq
  result$state <- out$state
  result$qvalues <- qvalues
  result$qvalues.rep <- qvalues.rep
  if(model.full == "statistical") {
    result$p.cutoff <- p.cutoff
  } else {
    result$pred1 <- pred1
    result$pred2 <- pred2
  }
  class(result) <- "convaq"
  
  return(result)
}
