#include "abspropagator.h"

const int ResPropType::NR       = 0;
const int ResPropType::RBW      = 1;
const int ResPropType::GS       = 2;
const int ResPropType::RhoOmega = 3;
const int ResPropType::Bugg     = 4;
const int ResPropType::VDst     = 5;
const int ResPropType::Flatte   = 6;

//AbsPropagator::AbsPropagator(const double& m, const double& p0, AbsVarWidth* width):
//  m_m(m),m_p0(p0),m_width(width)
//{
//}

AbsPropagator::AbsPropagator(const double& m, const double& p0):
  m_m(m),m_p0(p0)
{
}
