#include <algorithm>
#include "defines.h"
#include "Predicate.h"

bool comp_less(double a, double b) { return a < b; }
bool comp_greater(double a, double b) { return a > b; }
bool comp_leq(double a, double b) { return a <= b; }
bool comp_geq(double a, double b) { return a >= b; }

bool Predicate::match(const std::vector<std::vector<bool>> &state) {
  double freq;
  if(type == Normal) {
    int count = 0;
    for(size_t i = 0; i < state[0].size(); ++i) {
      bool any = false;
      for(size_t type = 0; type < 3; ++type) {
        any = any || state[type][i];
      }
      if(!any) ++count;
    }
    freq = (double)count / state[0].size();
  } else {
    freq = (double)std::accumulate(state[type].begin(), state[type].end(), 0) / state[type].size();
  }
  if(eq == EQ_NEQ) freq = 1.0 - freq;
  return f_comp(freq, value);
}

Predicate make_predicate(COMPARISON comp, double value, EQUALITY eq, VARIATION_TYPE type) {
  std::function<bool(double,double)> f_comp;
  if(comp == COMP_LESS) f_comp = comp_less;
  else if(comp == COMP_GREATER) f_comp = comp_greater;
  else if(comp == COMP_LEQ) f_comp = comp_leq;
  else if(comp == COMP_GEQ) f_comp = comp_geq;
  
  return Predicate(f_comp, value, eq, type);
}