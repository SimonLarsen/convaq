#include <Rcpp.h>
#include <vector>
#include "Segment.h"

using namespace Rcpp;

void df_to_segments(DataFrame df, std::vector<Segment> &segments) {
  IntegerVector patient = df["patient"];
  StringVector chr = df["chr"];
  IntegerVector start = df["start"];
  IntegerVector end = df["end"];
  IntegerVector type = df["type"];

  for(size_t i = 0; i < df.nrows(); ++i) {
    std::string chr_str(chr[i]);
    segments.emplace_back(patient[i], chr_str, start[i], end[i], type[i]);
  }
}
