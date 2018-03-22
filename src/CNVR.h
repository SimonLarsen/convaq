#ifndef CNVR_H
#define CNVR_H

#include "defines.h"
#include <string>

class CNVR {
public:
  std::string chr;
  int start;
  int end;
  int length;
  int type;
  double pvalue;
  double qvalue;

  CNVR(
    const std::string &chr,
    int start,
    int end,
    int length,
    int type,
    double pvalue,
    double qvalue = 0
  )
    : chr(chr),
      start(start),
      end(end),
      length(length),
      type(type),
      pvalue(pvalue),
      qvalue(qvalue)
  { }
};

#endif
