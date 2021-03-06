/** Copyright 2017 Vitaly Vorobyev
 ** @file buggwidth.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/buggwidth.h"

#include <cmath>

using std::sqrt;
using std::exp;
using std::pow;
using std::fabs;

const double BuggWidth::m_mr = 0.953;
const double BuggWidth::m_mrsq = pow(m_mr, 2);
const double BuggWidth::m_PI_Mass = 0.13957018;
const double BuggWidth::m_K_Mass = 0.493667;
const double BuggWidth::m_eta_Mass = 0.547862;
const double BuggWidth::m_sA = 0.41*pow(m_PI_Mass, 2);
const double BuggWidth::m_b1 = 1.302;
const double BuggWidth::m_b2 = 0.340;
const double BuggWidth::m_A = 2.426;
const double BuggWidth::m_g4pi = 0.011;
const double BuggWidth::m_alpha = 1.3;

BuggWidth::BuggWidth(void) : AbsVarWidth(1, 1, 1) {}

double BuggWidth::operator()(const double& s, const double& p) const {
    return s+p;
}

double BuggWidth::rho_4pi(const double& s) {
    return 1./(1.+exp(7.082-2.845*s));
}

double BuggWidth::rho(const double& m, const double& s) {
    const double var = 1.-4.*m*m/s;
    return var > 0 ? sqrt(var) : 0;
}

double BuggWidth::j1(const double& s) {
    return 1./M_PI*(2.+m_rho1_ps*log((1-m_rho1_ps)/(1+m_rho1_ps)));
}

double BuggWidth::z(const double& s) {
  return j1(s) - j1(m_mrsq);
}

double BuggWidth::g1sq(const double& s) {
  return m_mr*(m_b1+m_b2*s)*exp(-(s-m_mrsq)/m_A);
}

double BuggWidth::mrGamma1(const double& s) {
  return m_g1sq_pc*(s-m_sA)/(m_mrsq-m_sA)*m_rho1_ps;
}

double BuggWidth::mrGamma2(const double& s) {
  return 0.6*m_g1sq_pc*s/m_mrsq*
          exp(-m_alpha*fabs(s-4.*m_K_Mass*m_K_Mass))*rho(m_K_Mass, s);
}

double BuggWidth::mrGamma3(const double& s) {
  return 0.2*m_g1sq_pc*s/m_mrsq*
          exp(-m_alpha*fabs(s-4.*m_eta_Mass*m_eta_Mass))*rho(m_eta_Mass, s);
}

double BuggWidth::mrGamma4(const double& s) {
    return m_mr*m_g4pi*rho_4pi(s)/rho_4pi(m_mrsq);
}

void BuggWidth::GetWidths(const double& s, double* G1, double* GTot) {
    m_g1sq_pc = g1sq(s);
    m_rho1_ps = rho(m_PI_Mass, s);
    m_z_pc = z(s);
    *G1 = mrGamma1(s);
    *GTot = *G1 + mrGamma2(s) + mrGamma3(s) + mrGamma4(s);
}
