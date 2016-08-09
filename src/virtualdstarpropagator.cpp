#include "virtualdstarpropagator.h"
#include "constwidth.h"

using namespace std;

VirtualDstarPropagator::VirtualDstarPropagator(const double &beta1, const double &beta2):
  AbsPropagator(2.01,0),m_b1(beta1),m_b2(beta2)
{
}

compld VirtualDstarPropagator::operator()(const double& s,const double& p) const{
  return exp(-(m_b1+imone*m_b2)*s);
}
