#ifndef REGION_H
#define REGION_H

#include <string>
#include <vector>
#include <algorithm>

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
  
  double get_freq(size_t group, uint type) {
    int count = std::accumulate(state[group][type].begin(), state[group][type].end(), 0);
    return (double)count / state[0][0].size();
  }
};

#endif
