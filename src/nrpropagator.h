#ifndef NRPROPAGATOR_H
#define NRPROPAGATOR_H

#include "abspropagator.h"

class NRPropagator : public AbsPropagator{
public:
  NRPropagator(const double& alpha);
  EvtComplex operator()(const double& s,const double& p=0) const;

private:
  double m_alpha;
};

#endif // NRPROPAGATOR_H
