/** Copyright 2017 Vitaly Vorobyev
 ** @file virtualdstarpropagator.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include <cmath>

#include "../src/virtualdstarpropagator.h"
#include "../src/constwidth.h"

const std::complex<double> imone(0, 1);

VirtualDstarPropagator::VirtualDstarPropagator(const double &beta1,
                                               const double &beta2) :
    AbsPropagator(2.01, 0), m_b1(beta1), m_b2(beta2) {}

std::complex<double> VirtualDstarPropagator::operator()(
        const double& s, const double& p) const {
    return std::exp(-(m_b1 + imone * m_b2) * s);
}
