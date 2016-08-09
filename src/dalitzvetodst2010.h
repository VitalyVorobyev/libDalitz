#ifndef DALITZVETODST2010_H
#define DALITZVETODST2010_H

#include "dalitzveto.h"

class DalitzVetoDst2010 : public DalitzVeto{
public:
  DalitzVetoDst2010(const double cut);
  bool operator()(const double& AB,const double& BC = 0) const;

private:
  double m_mdst;
  double m_cut;
};

#endif // DALITZVETODST2010_H
