#ifndef VIRTUALDSTARPROPAGATOR2_H
#define VIRTUALDSTARPROPAGATOR2_H

#include "abspropagator.h"
#include "consts.h"

class VirtualDstarPropagator2 : public AbsPropagator{
public:
  VirtualDstarPropagator2();
  compld operator()(const double& s, const double& p) const;

private:
  double m_r;
};

#endif // VIRTUALDSTARPROPAGATOR2_H
