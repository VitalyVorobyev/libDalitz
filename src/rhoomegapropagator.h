#ifndef RHOOMEGAPROPAGATOR_H
#define RHOOMEGAPROPAGATOR_H

#include "gounarissakurai.h"
#include "relbreitwigner.h"
#include "consts.h"

/// \brief Class implementing rho(770)-omega(782) interfering propagator.
/// Taken from Eq.(13) in Phys. Rev. D92, 032002 (2015).
/// Introduced for the first time in R. R. Akhmetshin et al. (CMD-2 Collaboration),
/// Measurement of e+e- -> pi+ pi- cross section with CMD-2 around rho-meson,
/// Phys. Lett. B 527, 161 (2002).

class RhoOmegaPropagator : public GounarisSakurai{
public:
  RhoOmegaPropagator(const double& a, const double& theta);

  compld operator()(const double& s, const double& p) const;

private:
  compld m_omega_amp;
  RelBreitWigner* omega_prop;
};

#endif // RHOOMEGAPROPAGATOR_H
