#include "gswidth.h"
#include "math.h"

GSWidth::GSWidth(const double &G0, const double &m, const double &p0):
  AbsVarWidth(G0,m,p0), m_precalc(G0*m/pow(p0,3))
{
}

double GSWidth::operator()(const double& s, const double& p) const{
  return m_precalc*pow(p,3)/sqrt(s);
}
