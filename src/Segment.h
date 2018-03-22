#ifndef SEGMENT_H
#define SEGMENT_H

#include <string>

class Segment {
public:
  int patient;
  std::string chr;
  int start;
  int end;
  int type;

  Segment(int patient, const std::string &chr, int start, int end, int type)
    : patient(patient),
      chr(chr),
      start(start),
      end(end),
      type(type)
  {}
};

#endif
