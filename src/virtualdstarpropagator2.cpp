#include "virtualdstarpropagator2.h"
#include "blattweisskopf.h"

using namespace std;

VirtualDstarPropagator2::VirtualDstarPropagator2():
  AbsPropagator(2.01,0), m_r(BlattWeisskopf::m_r_meson)
{
}

compld VirtualDstarPropagator2::operator()(const double& s,const double& p) const{
  return imone;//std::exp(-(m_b1+ione*m_b2)*s);
}
