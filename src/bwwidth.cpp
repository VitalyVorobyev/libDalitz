#include "bwwidth.h"
#include "math.h"

BWWidth::BWWidth(const double& G0, const double &m, const double &p0, const int mom):
  AbsVarWidth(G0,m,p0),m_mom(mom)
{
  m_precalc = G0*m/pow(p0,2*mom+1);
}

double BWWidth::operator()(const double& s, const double& p) const{
  const double ff = (*m_ff)(p);
  return m_precalc*ff*ff*pow(p,2*m_mom+1)/sqrt(s);
}
