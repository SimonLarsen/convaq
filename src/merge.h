#ifndef MERGE_H
#define MERGE_H

#include <vector>
#include "CNVR.h"

void merge_adjacent(std::vector<CNVR> &regions, unsigned int threshold);

#endif
