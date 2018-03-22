#ifndef CNVR_H
#define CNVR_H

#include "defines.h"
#include <string>

class CNVR {
public:
  uint region;
  std::string chr;
  int start;
  int end;
  int length;
  int type;
  double pvalue;
  double qvalue;

  CNVR(
    uint region,
    const std::string &chr,
    int start,
    int end,
    int length,
    int type,
    double pvalue,
    double qvalue = 0
  )
    : region(region),
      chr(chr),
      start(start),
      end(end),
      length(length),
      type(type),
      pvalue(pvalue),
      qvalue(qvalue)
  { }
};

#endif
