#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "Region.h"
#include "CNVR.h"
#include "fisher_test.h"

void statistical_model(const std::vector<Region> &regions, int npatients1, int npatients2, double cutoff, std::vector<CNVR> &result) {
  for(size_t i = 0; i < regions.size(); ++i) {
    const Region &r = regions[i];
    for(size_t type = 0; type < 3; ++type) {
      int pos1 = std::accumulate(r.state[0][type].begin(), r.state[0][type].end(), 0);
      int neg1 = npatients1 - pos1;
      int pos2 = std::accumulate(r.state[1][type].begin(), r.state[1][type].end(), 0);
      int neg2 = npatients2 - pos2;

      double pval = fisher_test(pos1, neg1, pos2, neg2);
      if(pval <= cutoff) {
        result.emplace_back(i, r.chr, r.start, r.end, r.length, type, pval);
      }
    }
  }
}
