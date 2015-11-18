#include "flattewidth.h"
#include "consts.h"
#include "math.h"

FlatteWidth::FlatteWidth(const double &m):
  AbsVarWidth(1,m,1),m_g1(199),m_g2(m_g1*3.0),
  m_pi_sq(m_PI_Mass*m_PI_Mass),
  m_pi0_sq(m_PI0_Mass*m_PI0_Mass),
  m_K_sq(m_K_Mass*m_K_Mass),
  m_K0_sq(m_Ks0_Mass*m_Ks0_Mass)
{
}

double FlatteWidth::operator()(const double& s, const double& p) const{
  return m_g1*rho_pipi(s) + m_g2*rho_KK(s);
}

double FlatteWidth::rho_pipi(const double& s) const{
  return 2./3.*phsp_factor(m_pi_sq,s)+1./3.*phsp_factor(m_pi0_sq,s);
}
double FlatteWidth::rho_KK(const double& s) const{
  return 0.5*phsp_factor(m_K_sq,s)+0.5*phsp_factor(m_K0_sq,s);
}

double FlatteWidth::phsp_factor(const double& msq, const double&s) const{
  return sqrt(1.-4.*msq/s);
}
