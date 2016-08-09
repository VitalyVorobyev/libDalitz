#ifndef NRPROPAGATOR_H
#define NRPROPAGATOR_H

#include "abspropagator.h"
#include "consts.h"

class NRPropagator : public AbsPropagator{
public:
  NRPropagator(const double& alpha);
  compld operator()(const double& s,const double& p=0) const;

private:
  double m_alpha;
};

#endif // NRPROPAGATOR_H
