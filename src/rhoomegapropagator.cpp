/** Copyright 2017 Vitaly Vorobyev
 ** @file rhoomegapropagator.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/rhoomegapropagator.h"

#include <cmath>

using std::sqrt;
using std::pow;

const double RhoOmegaPropagator::m_PI_Mass = 0.13957018;
const double RhoOmegaPropagator::m_rho770_Mass = 0.77526;
const double RhoOmegaPropagator::m_rho770_Width = 0.1478;
const double RhoOmegaPropagator::m_omega_Mass = 0.78265;
const double RhoOmegaPropagator::m_omega_Width = 0.00849;
const double RhoOmegaPropagator::m_p0_rho770 =
        sqrt(0.25 * pow(pow(m_rho770_Mass, 2) - 2.*pow(m_PI_Mass, 2), 2) -
                    pow(pow(m_PI_Mass, 2) / m_rho770_Mass, 2));
const double RhoOmegaPropagator::m_p0_omega =
        sqrt(0.25 * pow(pow(m_omega_Mass, 2) - 2.*pow(m_PI_Mass, 2), 2) -
                    pow(pow(m_PI_Mass, 2) / m_omega_Mass, 2));

RhoOmegaPropagator::RhoOmegaPropagator(const double &a, const double &theta) :
    GounarisSakurai(m_rho770_Width, m_rho770_Mass, m_p0_rho770),
    m_omega_amp(a*std::complex<double>(cos(theta), sin(theta))) {
    omega_prop = new RelBreitWigner(m_omega_Width, m_omega_Mass, m_p0_omega,
                                    1, VarWType::Const);
}

std::complex<double> RhoOmegaPropagator::operator()(const double& s,
                                                    const double& p) const {
    return GounarisSakurai::operator ()(s, p) *
            (1. + m_omega_amp * (*omega_prop)(s, p));
}
