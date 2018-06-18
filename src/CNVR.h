#ifndef CNVR_H
#define CNVR_H

#include "defines.h"
#include "Region.h"
#include <string>
#include <numeric>
#include <set>

class CNVR {
public:
  std::string chr;
  int start;
  int end;
  int length;
  int type;
  double pvalue;
  double qvalue;
  std::vector<const Region*> regions;

  CNVR(const std::vector<CNVR> &_regions) {
    chr = _regions[0].chr;
    start = _regions[0].start;
    end = _regions[0].end;
    type = _regions[0].type;
    pvalue = _regions[0].pvalue;
    qvalue = _regions[0].qvalue;
    
    for(const CNVR &r : _regions) {
      start = std::min(start, r.start);
      end = std::max(end, r.end);
      pvalue = std::max(pvalue, r.pvalue);
      qvalue = std::max(qvalue, r.qvalue);
      for(const Region *re : r.regions) {
        regions.push_back(re);
      }
    }
    length = end-start+1;
  }
  
  CNVR(
    const Region &region,
    int type,
    double pvalue,
    double qvalue = 0
  )
    : chr(region.chr),
      start(region.start),
      end(region.end),
      length(region.length),
      type(type),
      pvalue(pvalue),
      qvalue(qvalue)
  {
    regions.push_back(&region);
  }
  
  std::vector<double> get_freq(size_t group, size_t type) const {
    std::vector<double> freq;
    for(const Region *r : regions) {
      int sum = std::accumulate(r->state[group][type].begin(), r->state[group][type].end(), 0);
      double f = (double)sum / r->state[group][type].size();
      freq.push_back(f);
    }
    std::sort(freq.begin(), freq.end());
    return(freq);
  }
  
  std::vector<std::vector<unsigned int>> get_state(size_t group) {
    size_t npatients = regions[0]->state[group][0].size();
    std::vector<std::vector<unsigned int>> state;
    for(size_t i = 0; i < npatients; ++i) {
      std::set<unsigned int> found;
      for(const Region *r : regions) {
        bool is_normal = true;
        for(size_t type = 0; type < 3; ++type) {
          if(r->state[group][type][i]) {
            found.insert(type);
            is_normal = false;
          }
        }
        if(is_normal) found.insert(Normal);
      }
      state.emplace_back(found.begin(), found.end());
    }
    return state;
  }
};

#endif
