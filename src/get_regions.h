#ifndef GET_REGIONS_H
#define GET_REGIONS_H

#include <string>
#include <vector>
#include <unordered_set>
#include "Segment.h"
#include "Region.h"
#include "Event.h"

void get_events(const std::vector<Segment> &segments, const std::string &chr, int group, std::vector<Event> &events);

void get_regions_chr(
    const std::vector<Segment> &segments1, const std::vector<Segment> &segments2,
    int npatients1, int npatients2,
    const std::string &chr,
    std::vector<Region> &regions
);

void get_regions(
  const std::vector<Segment> &segments1,
  const std::vector<Segment> &segments2,
  int npatients1,
  int npatients2,
  const std::unordered_set<std::string> &chromosomes,
  std::vector<Region> &regions
);

#endif
