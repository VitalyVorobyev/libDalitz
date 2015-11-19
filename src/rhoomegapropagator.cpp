#include "rhoomegapropagator.h"
#include "consts.h"

#include "math.h"

const double m_p0_rho770 = sqrt(0.25*pow(m_rho770_Mass*m_rho770_Mass-2.*m_PI_Mass*m_PI_Mass,2) - pow(m_PI_Mass*m_PI_Mass/m_rho770_Mass,2));
const double m_p0_omega  = sqrt(0.25*pow(m_omega_Mass*m_omega_Mass-2.*m_PI_Mass*m_PI_Mass,2) - pow(m_PI_Mass*m_PI_Mass/m_omega_Mass,2));

RhoOmegaPropagator::RhoOmegaPropagator(const double &a, const double &theta):
  GounarisSakurai(m_rho770_Width,m_rho770_Mass,m_p0_rho770),
  m_omega_amp(a*EvtComplex(cos(theta),sin(theta)))
{
  omega_prop = new RelBreitWigner(m_omega_Width,m_omega_Mass,m_p0_omega,1,VarWType::Const);
}

EvtComplex RhoOmegaPropagator::operator()(const double& s, const double& p) const{
  return GounarisSakurai::operator ()(s,p) * (1.+m_omega_amp*(*omega_prop)(s,p));
}
