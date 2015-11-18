#ifndef RELBREITWIGNER_H
#define RELBREITWIGNER_H

#include "abspropagator.h"
#include "bwwidth.h"
#include "constwidth.h"

class RelBreitWigner : public AbsPropagator{
public:
  RelBreitWigner(const double& G0, const double& m, const double& p0, const int mom, const bool constwidth = false);
  ~RelBreitWigner();

  EvtComplex operator()(const double& s, const double& p) const;
private:
  const bool m_const_width;
};

#endif // RELBREITWIGNER_H
