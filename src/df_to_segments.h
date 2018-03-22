#ifndef DF_TO_SEGMENTS_H
#define DF_TO_SEGMENTS_H

#include <Rcpp.h>
#include <vector>
#include "Segment.h"

void df_to_segments(Rcpp::DataFrame df, std::vector<Segment> &segments);

#endif
