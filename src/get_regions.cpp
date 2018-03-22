#include <vector>
#include <unordered_set>
#include <algorithm>

#include "get_regions.h"
#include "Event.h"
#include "Segment.h"

void get_events(const std::vector<Segment> &segments, const std::string &chr, int group, std::vector<Event> &events) {
  for(const Segment &s : segments) {
    if(s.chr == chr) {
      events.emplace_back(s.patient, s.chr, s.start, s.type, true, group);
      events.emplace_back(s.patient, s.chr, s.end, s.type, false, group);
    }
  }
}

void get_regions_chr(
    const std::vector<Segment> &segments1, const std::vector<Segment> &segments2,
    int npatients1, int npatients2,
    const std::string &chr,
    std::vector<Region> &regions
) {
  std::vector<Event> events;

  get_events(segments1, chr, 0, events);
  get_events(segments2, chr, 1, events);

  std::sort(events.begin(), events.end(), [](Event &a, Event &b) { return a.position < b.position; });

  std::vector<std::vector<std::vector<bool>>> states(2);
  for(size_t i = 0; i < 2; ++i) {
    states[i].resize(3);
    for(size_t j = 0; j < 3; ++j) {
      if(i == 0) states[i][j].resize(npatients1, false);
      else states[i][j].resize(npatients2, false);
    }
  }

  int currentPos = 0;
  int nextPos = events[0].position;
  int i = 0;
  while(i < events.size()) {
    while(i < events.size() && events[i].position == nextPos) {
      Event &e = events[i];
      states[e.group][e.type][e.patient] = e.isStart;
      ++i;
    }
    if(i >= events.size()) break;

    currentPos = nextPos;
    nextPos = events[i].position;

    regions.emplace_back(chr, currentPos, nextPos, nextPos-currentPos+1, states);
  }
}

void get_regions(
  const std::vector<Segment> &segments1,
  const std::vector<Segment> &segments2,
  int npatients1,
  int npatients2,
  const std::unordered_set<std::string> &chromosomes,
  std::vector<Region> &regions
) {
  regions.clear();
  for(const std::string &chr : chromosomes) {
    get_regions_chr(segments1, segments2, npatients1, npatients2, chr, regions);
  }
}
