#ifndef RHOOMEGAPROPAGATOR_H
#define RHOOMEGAPROPAGATOR_H

#include "gounarissakurai.h"
#include "relbreitwigner.h"

class RhoOmegaPropagator : public GounarisSakurai{
public:
  RhoOmegaPropagator(const double& a, const double& theta);

  EvtComplex operator()(const double& s, const double& p) const;

private:
  EvtComplex m_omega_amp;
  RelBreitWigner* omega_prop;
};

#endif // RHOOMEGAPROPAGATOR_H
