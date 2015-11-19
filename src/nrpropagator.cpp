#include "nrpropagator.h"
#include "constwidth.h"

#include <complex>

NRPropagator::NRPropagator(const double &alpha):
  AbsPropagator(2.01,0,new ConstWidth(0.1)),
  m_alpha(alpha)
{
}

EvtComplex NRPropagator::operator()(const double& s,const double& p) const{
  const std::complex<double> ione(0, 1);
  const std::complex<double> res = std::exp(ione*m_alpha*s);
  return EvtComplex(res.real(),res.imag());
}
