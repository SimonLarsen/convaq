#ifndef EVENT_H
#define EVENT_H

#include <string>

class Event {
public:
  int patient;
  std::string chr;
  int position;
  int type;
  bool isStart;
  int group;

  Event(int patient, const std::string &chr, int position, int type, bool isStart, int group)
    : patient(patient),
      chr(chr),
      position(position),
      type(type),
      isStart(isStart),
      group(group)
  {}

};

#endif
