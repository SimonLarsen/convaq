#ifndef STATISTICAL_MODEL
#define STATISTICAL_MODEL

#include <vector>
#include "CNVR.h"
#include "Region.h"

void statistical_model(const std::vector<Region> &regions, int npatients1, int npatients2, double cutoff, std::vector<CNVR> &result);

#endif
