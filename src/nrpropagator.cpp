/** Copyright 2017 Vitaly Vorobyev
 ** @file nrpropagator.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/nrpropagator.h"

#include <cmath>

NRPropagator::NRPropagator(const double &alpha) :
  AbsPropagator(2.01, 1), m_alpha(alpha) {}

std::complex<double> NRPropagator::operator()(const double& s,
                                              const double& p) const {
  if (m_alpha == 0) return 1.;
  return std::exp(std::complex<double>(0, 1) * m_alpha * s);
}
