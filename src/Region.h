#ifndef REGION_H
#define REGION_H

#include <string>
#include <vector>

class Region {
public:
  std::string chr;
  int start;
  int end;
  int length;
  std::vector<std::vector<std::vector<bool>>> state;

  Region(const std::string &chr, int start, int end, int length, const std::vector<std::vector<std::vector<bool>>> &state)
    : chr(chr),
      start(start),
      end(end),
      length(length),
      state(state)
  {}
};

#endif
