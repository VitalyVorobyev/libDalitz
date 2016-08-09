#ifndef BUGGPROPAGATOR_H
#define BUGGPROPAGATOR_H

#include "abspropagator.h"
#include "buggwidth.h"
#include "consts.h"

class BuggPropagator : public AbsPropagator{
public:
  BuggPropagator(void);
  ~BuggPropagator(void);

  compld operator()(const double& s, const double& p = 0) const;

private:
  BuggWidth* m_width;
};

#endif // BUGGPROPAGATOR_H
