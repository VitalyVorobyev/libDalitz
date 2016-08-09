#include "nrpropagator.h"
#include "constwidth.h"

using namespace std;

NRPropagator::NRPropagator(const double &alpha):
  AbsPropagator(2.01,1),m_alpha(alpha)
{
}

compld NRPropagator::operator()(const double& s,const double& p) const{
  if(m_alpha == 0) return 1.;
//  std::cout << "NRPropagator::operator()" << std::endl;
  return std::exp(imone*m_alpha*s);
}
