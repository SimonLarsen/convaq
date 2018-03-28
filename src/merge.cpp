#include <Rcpp.h>
#include <vector>
#include "merge.h"
#include "CNVR.h"

void merge_adjacent(std::vector<CNVR> &regions, uint threshold) {
  std::sort(regions.begin(), regions.end(), [](const CNVR &a, const CNVR &b){
    if(a.type != b.type) return a.type < b.type;
    if(a.chr != b.chr) return a.chr < b.chr;
    return(a.start < b.start);
  });
  
  size_t i = 0;
  std::vector<CNVR> out;
  while(i < regions.size()) {
    size_t first = i;
    size_t last = i++;
    
    while(
      i < regions.size() &&
      regions[i].type == regions[first].type &&
      regions[i].chr == regions[first].chr &&
      regions[i].start-regions[i-1].end-1 <= (int)threshold
    ) {
      last = i;
      ++i;
    }
    
    if(first == last) {
      out.push_back(regions[first]);
    } else {
      std::vector<CNVR> merge_regions(regions.begin()+first, regions.begin()+last+1);
      out.emplace_back(merge_regions);
    }
  }
  
  regions.clear();
  std::copy(out.begin(), out.end(), std::back_inserter(regions));
}