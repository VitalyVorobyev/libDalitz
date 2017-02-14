/** Copyright 2017 Vitaly Vorobyev
 ** @file buggproparator.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include <cmath>

#include "../src/buggpropagator.h"
#include "./buggwidth.h"

BuggPropagator::BuggPropagator(void) :
    AbsPropagator(0, 0), m_width(new BuggWidth()) {}

std::complex<double> BuggPropagator::operator()(const double& s,
                                                const double& p) const {
  double G1, GTot;
  m_width->GetWidths(s, &G1, &GTot);
  const double sA   = m_width->sA();
  const double mrsq = m_width->mrsq();
  const double z    = m_width->z();
  const double g1sq = m_width->g1sq();

  return G1/(mrsq-s-g1sq*(s-sA)*z/(mrsq-sA)-
             std::complex<double>(0, 1)*std::sqrt(mrsq)*GTot);
}

BuggPropagator::~BuggPropagator(void) {
    delete m_width;
}
