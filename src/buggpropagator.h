#ifndef BUGGPROPAGATOR_H
#define BUGGPROPAGATOR_H

#include "abspropagator.h"

class BuggPropagator : public AbsPropagator{
public:
  BuggPropagator(void);

  EvtComplex operator()(const double& s, const double& p = 0) const;
};

#endif // BUGGPROPAGATOR_H
