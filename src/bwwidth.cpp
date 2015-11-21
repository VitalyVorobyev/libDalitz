#include "bwwidth.h"
#include "math.h"

#include <iostream>

BWWidth::BWWidth(const double& G0, const double &m, const double &p0, const int mom):
  AbsVarWidth(G0,m,p0),m_mom(mom)
{
  m_precalc = G0*m/pow(p0,2*mom+1);
//  std::cout << "BWWidth1: m = " << m << " " << p0 << " " << m_precalc << std::endl;
  m_ff = new BlattWeisskopf(mom,p0,FFType::FFResonance);
}

double BWWidth::operator()(const double& s, const double& p) const{
  const double ff = (*m_ff)(p);
//  std::cout << "BWWidth2: s = " << s << " " << p << std::endl;
//  std::cout << "BWWidth::operator(): " << ff << " " << m_precalc*ff*ff*pow(p,2*m_mom+1)/sqrt(s);
//  std::cout << " " << G0() << std::endl;
  return m_precalc*ff*ff*pow(p,2*m_mom+1)/sqrt(s);
}
