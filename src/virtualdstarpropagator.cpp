#include "virtualdstarpropagator.h"
#include "constwidth.h"

#include <complex>

VirtualDstarPropagator::VirtualDstarPropagator(const double &beta1, const double &beta2):
  AbsPropagator(2.01,0,new ConstWidth(0.1)),
  m_b1(beta1),m_b2(beta2)
{
}

VirtualDstarPropagator::~VirtualDstarPropagator(void){
  delete (ConstWidth*) width();
}

EvtComplex VirtualDstarPropagator::operator()(const double& s,const double& p) const{
  const std::complex<double> ione(0, 1);
  const std::complex<double> res = std::exp(-(m_b1+ione*m_b2)*s);
  return EvtComplex(res.real(),res.imag());
}
