/** Copyright 2017 Vitaly Vorobyev
 ** @file virtualdstarpropagator2.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/virtualdstarpropagator2.h"
#include "../src/blattweisskopf.h"

VirtualDstarPropagator2::VirtualDstarPropagator2() :
    AbsPropagator(2.01, 0), m_r(BlattWeisskopf::m_r_meson) {}

std::complex<double> VirtualDstarPropagator2::operator()(
        const double& s, const double& p) const {
    return std::complex<double>(0, 1);
}
