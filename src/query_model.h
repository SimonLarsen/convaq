#ifndef QUERY_MODEL_H
#define QUERY_MODEL_H

#include <vector>
#include "defines.h"
#include "Region.h"
#include "CNVR.h"

void query_model(
    const std::vector<Region> &regions,
    int npatients1, int npatients2,
    COMPARISON comp1, double value1, EQUALITY eq1, VARIATION_TYPE type1,
    COMPARISON comp2, double value2, EQUALITY eq2, VARIATION_TYPE type2,
    std::vector<CNVR> &result
);

#endif