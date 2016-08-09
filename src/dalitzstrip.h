#ifndef DALITZSTRIP_H
#define DALITZSTRIP_H

#include "dalitzveto.h"

class DalitzStrip : public DalitzVeto{
public:
  DalitzStrip();

private:
  int m_type;
};

#endif // DALITZSTRIP_H
