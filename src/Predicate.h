#ifndef PREDICATE_H
#define PREDICATE_H

#include <vector>
#include <functional>
#include "defines.h"

class Predicate {
public:
  std::function<bool(double,double)> f_comp;
  double value;
  EQUALITY eq;
  VARIATION_TYPE type; 
  
  Predicate(
    std::function<bool(double,double)> f_comp,
    double value,
    EQUALITY eq,
    VARIATION_TYPE type
  )
    : f_comp(f_comp),
      value(value),
      eq(eq),
      type(type)
  {}
    
  bool match(const std::vector<std::vector<bool>> &state);
};

Predicate make_predicate(COMPARISON comp, double value, EQUALITY eq, VARIATION_TYPE type);

#endif