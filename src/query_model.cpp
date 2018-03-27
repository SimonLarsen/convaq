#include <Rcpp.h>
#include <vector>
#include "query_model.h"
#include "Region.h"
#include "CNVR.h"
#include "Predicate.h"

void query_model(
    const std::vector<Region> &regions,
    int npatients1, int npatients2,
    COMPARISON comp1, double value1, EQUALITY eq1, VARIATION_TYPE type1,
    COMPARISON comp2, double value2, EQUALITY eq2, VARIATION_TYPE type2,
    std::vector<CNVR> &result
) {
  Predicate pred1 = make_predicate(comp1, value1, eq1, type1);
  Predicate pred2 = make_predicate(comp2, value2, eq2, type2);
  
  for(size_t i = 0; i < regions.size(); ++i) {
    const Region &r = regions[i];
    if(pred1.match(r.state[0]) && pred2.match(r.state[1])) {
      result.emplace_back(r, Normal, 1);
    }
  }
}