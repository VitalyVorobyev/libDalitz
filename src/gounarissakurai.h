#ifndef GOUNARISSAKURAI_H
#define GOUNARISSAKURAI_H

#include "abspropagator.h"
#include "gswidth.h"
#include "constwidth.h"

/// \brief Class implementing Gounaris-Sakurai lineshape
/// Gounaris-Sakurai lineshape is introduced by
/// G. J. Gounaris and J. J. Sakurai in ''Finite Width Corrections to
/// the Vector-Meson-Dominance Prediction for rho -> e+e-``,
/// Phys. Rev. Lett. 21, 244 (1968).
///
/// GS lineshape is applied in LHCb experiment
/// (Eq.(11,12) in Phys. Rev. D92, 032002 (2015))

class GounarisSakurai : public AbsPropagator{
public:
  GounarisSakurai(const double& G0, const double& m, const double& p0, const bool constwidth = false);
  ~GounarisSakurai();

  EvtComplex operator()(const double& s, const double& p) const;

private:
  double f(const double& s, const double& p) const;
  double h(const double& s, const double& p) const;
  int hhder(double& h,double& hder,const double& p) const;
  double g(void);
  double m_g;
  const bool m_const_width;
  AbsVarWidth* m_width;
};

#endif // GOUNARISSAKURAI_H
