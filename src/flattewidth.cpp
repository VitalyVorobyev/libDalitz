/** Copyright 2017 Vitaly Vorobyev
 ** @file flattewidth.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/flattewidth.h"

#include <cmath>

using std::pow;

const double FlatteWidth::m_g1 = 0.199;
const double FlatteWidth::m_g2 = FlatteWidth::m_g1*3.0;
const double FlatteWidth::m_pi = 0.13957018;
const double FlatteWidth::m_pi_sq = pow(FlatteWidth::m_pi, 2);
const double FlatteWidth::m_pi0 = 0.1349766;
const double FlatteWidth::m_pi0_sq = pow(FlatteWidth::m_pi0, 2);
const double FlatteWidth::m_K = 0.493677;
const double FlatteWidth::m_K_sq = pow(FlatteWidth::m_K, 2);
const double FlatteWidth::m_K0 = 0.497611;
const double FlatteWidth::m_K0_sq = pow(FlatteWidth::m_K0, 2);

FlatteWidth::FlatteWidth(const double &m) : AbsVarWidth(1, m, 1) {}

double FlatteWidth::operator()(const double& s, const double& p) const {
    return m_g1*rho_pipi(s) + m_g2*rho_KK(s);
}

double FlatteWidth::rho_pipi(const double& s) const {
    return 2./3.*phsp_factor(m_pi_sq, s)+1./3.*phsp_factor(m_pi0_sq, s);
}
double FlatteWidth::rho_KK(const double& s) const {
    return 0.5*phsp_factor(m_K_sq, s)+0.5*phsp_factor(m_K0_sq, s);
}

double FlatteWidth::phsp_factor(const double& msq, const double&s) const {
    const double var = 1.-4.*msq/s;
    return var > 0 ? std::sqrt(var) : 0;
}
