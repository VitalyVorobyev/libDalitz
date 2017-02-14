/** Copyright 2017 Vitaly Vorobyev
 ** @file resdecayangulardistribution.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/resdecayangulardistribution.h"

#include <cmath>

#include "../src/dalitzphasespace.h"

using std::pow;

typedef ResDecayAngularDistribution RDAD;
typedef DalitzPhaseSpace DPhSp;

bool RDAD::m_use_mRsq = true;

RDAD::ResDecayAngularDistribution(const int spin, const double& mmo,
                                  const double& mca, const double& mcb,
                                  const double& mcc, const double& mres) :
    m_spin(spin), m_mmo(mmo), m_mca(mca), m_mcb(mcb), m_mcc(mcc),
    m_mmo_sq(mmo*mmo), m_mca_sq(mca*mca), m_mcb_sq(mcb*mcb),
    m_mcc_sq(mcc*mcc), m_mre_sq(mres*mres) {
    if      (m_spin == 1) make_precalc_one();
    else if (m_spin == 2) make_precalc_all();
}

double RDAD::operator()(const double& mACsq, const double& mBCsq) const {
    if (!m_spin) return 1.;
    double pdf = 1;
    const double mABsq = DPhSp::GetmBCsq(m_mmo_sq, m_mca_sq, m_mcb_sq,
                                         m_mcc_sq, mACsq, mBCsq);
    double cos_hel = m_spin > 2 ?
                DPhSp::CosHelAB(m_mmo, m_mca, m_mcb, m_mcc, mABsq) : 1;
    double p       = m_spin > 2 ? DPhSp::p(m_mre_sq, m_mmo, m_mcc) : 0;
    double q       = m_spin > 2 ? DPhSp::q(m_mre_sq, m_mca, m_mcb) : 0;
    double mRsq = m_mre_sq;
    switch (m_spin) {
    case 1:
        if (!m_use_mRsq) mRsq = mABsq;
        pdf = mACsq - mBCsq + m_pc_var1 / mRsq;
        break;
    case 2:
        if (!m_use_mRsq) mRsq = mABsq;
        pdf = pow(mACsq - mBCsq + m_pc_var1 / mRsq, 2) -
                (1. / 3.) * (mABsq - m_pc_var2 + m_pc_var3 / mRsq) *
                (mABsq - m_pc_var4 + m_pc_var5 / mRsq);
        break;
    case 3:  // spin 3 resonance Zemach helicity factor
        pdf = -8.*(5.*pow(cos_hel, 3) - 3.*cos_hel) / 5.*pow(p*q, 3);
        break;
    case 4:  // spin 4 resonance Zemach helicity factor
        pdf = 16.*(35.*pow(cos_hel, 4) - 30.*pow(cos_hel, 2) + 3.) /
                35.*pow(p*q, 4);
        break;
    case 5:  // spin 5 resonance Zemach helicity factor
        pdf = -32.*(63.*pow(cos_hel, 5) - 70.*pow(cos_hel, 3) + 15.*cos_hel) /
                63.*pow(p*q, 5);
    }
    return pdf;
}

void RDAD::make_precalc_one(void) {
    m_pc_var1 = (m_mmo_sq - m_mcc_sq) * (m_mcb_sq - m_mca_sq);
}

void RDAD::make_precalc_all(void) {
    make_precalc_one();
    m_pc_var2 = 2.*(m_mmo_sq+m_mcc_sq);
    m_pc_var3 = pow(m_mmo_sq-m_mcc_sq, 2);
    m_pc_var4 = 2.*(m_mca_sq+m_mcb_sq);
    m_pc_var5 = pow(m_mca_sq-m_mcb_sq, 2);
}
