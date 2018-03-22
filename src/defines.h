#ifndef DEFINES_H
#define DEFINES_H

enum VARIATION_TYPE {
  Gain   = 0,
  Loss   = 1,
  LOH    = 2,
  Normal = 3
};

enum MODEL {
  MODEL_STAT  = 1,
  MODEL_QUERY = 2
};

enum COMPARISON {
  COMP_NONE    = 0,
  COMP_LESS    = 1,
  COMP_GREATER = 2,
  COMP_LEQ     = 3,
  COMP_GEQ     = 4
};

enum EQUALITY {
  EQ_NONE  = 0,
  EQ_EQUAL = 1,
  EQ_NEQ   = 2
};

#endif
