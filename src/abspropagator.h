#ifndef ABSPROPAGATOR_H
#define ABSPROPAGATOR_H

#include "absvarwidth.h"
#include "EvtComplex.h"

class AbsPropagator{
public:
  AbsPropagator(const double &m, const double &p0, AbsVarWidth* width);
  virtual ~AbsPropagator() = 0;
  virtual EvtComplex operator()(const double& s, const double& p) const = 0;

  double m(void)  const {return m_m;}
  double p0(void) const {return m_p0;}
  AbsVarWidth* width(void) const {return m_width;}

private:
  double m_m;
  double m_p0;
  AbsVarWidth* m_width;
};

#endif // ABSPROPAGATOR_H
