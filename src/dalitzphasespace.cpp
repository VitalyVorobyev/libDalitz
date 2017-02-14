/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzphasespace.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/dalitzphasespace.h"
#include <cmath>
#include <iostream>

typedef DalitzPhaseSpace DPhSp;

using std::cout;
using std::endl;
using std::pow;
using std::sqrt;
using linal::LVect;

DPhSp::DalitzPhaseSpace(const double& mM, const double& mA,
                                   const double& mB, const double& mC) :
    m_mM(mM), m_mA(mA), m_mB(mB), m_mC(mC),
    m_mMsq(m_mM*m_mM), m_mAsq(m_mA*m_mA),
    m_mBsq(m_mB*m_mB), m_mCsq(m_mC*m_mC),
    m_mABsq_min((m_mA+m_mB)*(m_mA+m_mB)),
    m_mABsq_max((m_mM-m_mC)*(m_mM-m_mC)),
    m_mACsq_min((m_mA+m_mC)*(m_mA+m_mC)),
    m_mACsq_max((m_mM-m_mB)*(m_mM-m_mB)),
    m_mBCsq_min((m_mB+m_mC)*(m_mB+m_mC)),
    m_mBCsq_max((m_mM-m_mA)*(m_mM-m_mA)),
    m_mass_sq_sum(m_mMsq+m_mAsq+m_mBsq+m_mCsq),
    m_area(0) {}

DalitzPhaseSpace::DalitzPhaseSpace(const DalitzPhaseSpace& phsp) :
    DalitzPhaseSpace(phsp.mM(), phsp.mA(), phsp.mB(), phsp.mC()) {}

double DPhSp::GetmBCsq(const double& mABsq, const double& mACsq) const {
    const double mBCsq = m_mass_sq_sum-mABsq-mACsq;
    if (mBCsq < 0) {
        cout << "DPhSp::GetmBCsq: " << m_mass_sq_sum << " " << mABsq << " "
             << mACsq << " " << mBCsq << endl;
        cout << "IsInPlot: " << IsInPlot(mABsq, mACsq) << endl;
    }
    return mBCsq;
}

/// Min value of AB squared invariant mass
double DPhSp::GetmABsqMin(const double& mA, const double& mB) {
    return pow(mA+mB, 2);
}

/// Max value of AB squared invariant mass
double DPhSp::GetmABsqMax(const double& mM, const double& mC) {
    return pow(mM-mC, 2);
}

// Particle energies in the resonance rest frame
double DPhSp::eA_AB(const double& mABsq) const {
    return eA_AB(mABsq, m_mAsq, m_mBsq);
}
double DPhSp::eB_AB(const double& mABsq) const {
    return eB_AB(mABsq, m_mAsq, m_mBsq);
}
double DPhSp::eC_AB(const double& mABsq) const {
    return eC_AB(mABsq, m_mMsq, m_mCsq);
}
double DPhSp::eM_AB(const double& mABsq) const {
    return eM_AB(mABsq, m_mMsq, m_mCsq);
}

double DPhSp::eA_AC(const double& mACsq) const {
    return eA_AC(mACsq, m_mAsq, m_mCsq);
}
double DPhSp::eB_AC(const double& mACsq) const {
    return eB_AC(mACsq, m_mMsq, m_mBsq);
}
double DPhSp::eC_AC(const double& mACsq) const {
    return eC_AC(mACsq, m_mAsq, m_mCsq);
}
double DPhSp::eM_AC(const double& mACsq) const {
    return eM_AC(mACsq, m_mMsq, m_mBsq);
}

double DPhSp::eA_BC(const double& mBCsq) const {
    return eA_BC(mBCsq, m_mMsq, m_mAsq);
}
double DPhSp::eB_BC(const double& mBCsq) const {
    return eB_BC(mBCsq, m_mBsq, m_mCsq);
}
double DPhSp::eC_BC(const double& mBCsq) const {
    return eC_BC(mBCsq, m_mBsq, m_mCsq);
}
double DPhSp::eM_BC(const double& mBCsq) const {
    return eM_BC(mBCsq, m_mMsq, m_mAsq);
}

// Particle momenta in the resonance rest frame
double DPhSp::pA_AB(const double& mABsq) const {
    return pA_AB(mABsq, m_mAsq, m_mBsq);
}
double DPhSp::pB_AB(const double& mABsq) const {
    return pB_AB(mABsq, m_mAsq, m_mBsq);
}
double DPhSp::pC_AB(const double& mABsq) const {
    return pC_AB(mABsq, m_mMsq, m_mCsq);
}
double DPhSp::pM_AB(const double& mABsq) const {
    return pM_AB(mABsq, m_mMsq, m_mCsq);
}

double DPhSp::pA_AC(const double& mACsq) const {
    return pA_AC(mACsq, m_mAsq, m_mCsq);
}
double DPhSp::pB_AC(const double& mACsq) const {
    return pB_AC(mACsq, m_mMsq, m_mBsq);
}
double DPhSp::pC_AC(const double& mACsq) const {
    return pC_AC(mACsq, m_mAsq, m_mCsq);
}
double DPhSp::pM_AC(const double& mACsq) const {
    return pM_AC(mACsq, m_mMsq, m_mBsq);
}

double DPhSp::pA_BC(const double& mBCsq) const {
    return pA_BC(mBCsq, m_mMsq, m_mAsq);
}
double DPhSp::pB_BC(const double& mBCsq) const {
    return pB_BC(mBCsq, m_mBsq, m_mCsq);
}
double DPhSp::pC_BC(const double& mBCsq) const {
    return pC_BC(mBCsq, m_mBsq, m_mCsq);
}
double DPhSp::pM_BC(const double& mBCsq) const {
    return pM_BC(mBCsq, m_mMsq, m_mAsq);
}

// Resonance energy in the M frame
double DPhSp::eAB_M(const double& mABsq) const {
    return eAB_M(m_mMsq, mABsq, m_mCsq);
}
double DPhSp::eAC_M(const double& mACsq) const {
    return eAC_M(m_mMsq, mACsq, m_mBsq);
}
double DPhSp::eBC_M(const double& mBCsq) const {
    return eBC_M(m_mMsq, mBCsq, m_mAsq);
}

// Resonance momentum in the M frame
double DPhSp::pAB_M(const double& mABsq) const {
    return pAB_M(m_mMsq, mABsq, m_mCsq);
}
double DPhSp::pAC_M(const double& mACsq) const {
    return pAC_M(m_mMsq, mACsq, m_mBsq);
}
double DPhSp::pBC_M(const double& mBCsq) const {
    return pBC_M(m_mMsq, mBCsq, m_mAsq);
}

// Helicity angles
double DPhSp::CosHelAB(const double& mABsq, const double& mACsq,
                       double* pq) const {
    return CosHelAB(m_mMsq, m_mAsq, m_mBsq, m_mCsq, mABsq, mACsq, pq);
}
double DPhSp::CosHelAC(const double& mACsq, const double& mBCsq,
                       double* pq) const {
    return CosHelAC(m_mMsq, m_mAsq, m_mBsq, m_mCsq, mACsq, mBCsq, pq);
}
double DPhSp::CosHelBC(const double& mBCsq, const double& mABsq,
                       double* pq) const {
    return CosHelBC(m_mMsq, m_mAsq, m_mBsq, m_mCsq, mBCsq, mABsq, pq);
}

void DPhSp::RangeCalc(const double& e1, const double& e2, const double& m1sq,
                      const double& m2sq, double* msqmin, double* msqmax) {
    const double p1sq = e1*e1-m1sq;
    const double p1   = p1sq > 0 ? sqrt(p1sq) : 0;
    const double p2sq = e2*e2-m2sq;
    const double p2   = p2sq > 0 ? sqrt(p2sq) : 0;
    const double esum = e1+e2;
    const double pdif = p1-p2;
    const double psum = p1+p2;
    *msqmax = (esum-pdif)*(esum+pdif);
    *msqmin = (esum-psum)*(esum+psum);
}

int DPhSp::mABsqRange_AC(const double& mACsq,
                         double* mABsqMin, double* mABsqMax) const {
    if (mACsq > m_mACsq_max || mACsq < m_mACsq_min) return -1;
    RangeCalc(eA_AC(mACsq), eB_AC(mACsq), m_mAsq, m_mBsq, mABsqMin, mABsqMax);
    return 0;
}

int DPhSp::mBCsqRange_AC(const double& mACsq,
                         double* mBCsqMin, double* mBCsqMax) const {
    if (mACsq > m_mACsq_max || mACsq < m_mACsq_min) return -1;
    RangeCalc(eC_AC(mACsq), eB_AC(mACsq), m_mCsq, m_mBsq, mBCsqMin, mBCsqMax);
    return 0;
}

int DPhSp::mBCsqRange_AB(const double& mABsq,
                         double* mBCsqMin, double* mBCsqMax) const {
    if (mABsq > m_mABsq_max || mABsq < m_mABsq_min) return -1;
    RangeCalc(eB_AB(mABsq), eC_AB(mABsq), m_mBsq, m_mCsq, mBCsqMin, mBCsqMax);
    return 0;
}

int DPhSp::mACsqRange_AB(const double& mABsq,
                         double* mACsqMin, double* mACsqMax) const {
    if (mABsq > m_mABsq_max || mABsq < m_mABsq_min) return -1;
    RangeCalc(eA_AB(mABsq), eC_AB(mABsq), m_mAsq, m_mCsq, mACsqMin, mACsqMax);
    return 0;
}

int DPhSp::mACsqRange_BC(const double& mBCsq,
                         double* mACsqMin, double* mACsqMax) const {
    if (mBCsq > m_mBCsq_max || mBCsq < m_mBCsq_min) return -1;
    RangeCalc(eC_BC(mBCsq), eA_BC(mBCsq), m_mCsq, m_mAsq, mACsqMin, mACsqMax);
    return 0;
}

int DPhSp::mABsqRange_BC(const double& mBCsq,
                         double* mABsqMin, double* mABsqMax) const {
    if (mBCsq > m_mBCsq_max || mBCsq < m_mBCsq_min) return -1;
    RangeCalc(eB_BC(mBCsq), eA_BC(mBCsq), m_mBsq, m_mAsq, mABsqMin, mABsqMax);
    return 0;
}

bool DPhSp::IsInPlot(const double& mABsq, const double& mACsq) const {
    double mABsqmin, mABsqmax;
    if (mABsqRange_AC(mACsq, &mABsqmin, &mABsqmax)) return false;
    if (mABsq < mABsqmin || mABsq > mABsqmax) return false;
    return true;
}

// * Static methods * //
double DPhSp::GetmBCsq(const double& mMsq, const double& mAsq,
                       const double& mBsq, const double& mCsq,
                       const double& mABsq, const double& mACsq) {
    return mMsq+mAsq+mBsq+mCsq-mABsq-mACsq;
}

double DPhSp::pRes(const double& mMsq, const double& mAsq,
                   const double& mBsq) {
    const double res = (0.25*pow(mMsq-mAsq-mBsq, 2)-mAsq*mBsq)/mMsq;
    return res > 0 ? sqrt(res) : 0;
}

double DPhSp::q(const double& mRsq, const double& mA, const double& mB) {
    const double var = (mRsq-pow(mA+mB, 2))*(mRsq-pow(mA-mB, 2));
    return var > 0 ? sqrt(var/(4.*mRsq)) : 0;
}

double DPhSp::p(const double& mRsq, const double& mM, const double& mC) {
    const double var = (pow(mM-mC, 2)-mRsq)*(pow(mM+mC, 2)-mRsq);
    return var > 0 ? sqrt(var/(4.*mRsq)) : 0;
}

double DPhSp::ysq(const double& mMsq, const double& mABsq,
                  const double& mCsq) {
// ??? to be checked
    return 0.5*(mMsq + mABsq - mCsq)/sqrt(mABsq*mMsq) - 1.;
}

double DPhSp::CosHelAB(const double& mMo, const double& mA,
                       const double& mB, const double& mC,
                       const double& mABsq) {
    const double mABmin = GetmABsqMin(mA, mB);
    const double mABmax = GetmABsqMax(mMo, mC);
    return (mABmax+mABmin-2.*mABsq)/(mABmax-mABmin);
}

double DPhSp::CosHelAB(const double& mMsq, const double& mAsq,
                       const double& mBsq, const double& mCsq,
                       const double& mABsq, const double& mACsq, double* pq) {
    const double eA = eA_AB(mABsq, mAsq, mBsq);
    const double pAsq = pow(eA, 2) - mAsq;
    if (pAsq < 0) {*pq = 0; return 0;}
    const double pA = sqrt(pAsq);
    const double eC = eC_AB(mABsq, mMsq, mCsq);
    const double pCsq = pow(eC, 2) - mCsq;
    if (pCsq < 0) {*pq = 0; return 0;}
    const double pC = sqrt(pCsq);
    *pq = pA*pC;
    return (mAsq + mCsq + 2.*eA*eC - mACsq)/(2.*pA*pC);
}

double DPhSp::CosHelAC(const double& mMsq, const double& mAsq,
                       const double& mBsq, const double& mCsq,
                       const double& mACsq, const double& mBCsq, double* pq) {
    const double eB = eB_AC(mACsq, mMsq, mBsq);
    const double pBsq = pow(eB, 2) - mBsq;
    if (pBsq < 0) {*pq = 0; return 0;}
    const double pB = sqrt(pBsq);
    const double eC = eC_AC(mACsq, mAsq, mCsq);
    const double pCsq = pow(eC, 2) - mCsq;
    if (pCsq < 0) {*pq = 0; return 0;}
    const double pC   = sqrt(pCsq);
    *pq = pB*pC;
    return (mBsq + mCsq + 2.*eB*eC - mBCsq)/(2.*pB*pC);
}

double DPhSp::CosHelBC(const double& mMsq, const double& mAsq,
                       const double& mBsq, const double& mCsq,
                       const double& mBCsq, const double& mABsq, double* pq) {
    const double eA = eA_BC(mBCsq, mMsq, mAsq);
    const double pAsq = pow(eA, 2) - mAsq;
    if (pAsq < 0) {*pq = 0; return 0;}
    const double pA = sqrt(pAsq);
    const double eB = eB_BC(mBCsq, mBsq, mCsq);
    const double pBsq = pow(eB, 2) - mBsq;
    if (pBsq < 0) {*pq = 0; return 0;}
    const double pB = sqrt(pBsq);
    *pq = pA*pB;
    return (mAsq + mBsq + 2.*eA*eB - mABsq)/(2.*pA*pB);
}

double DPhSp::GetEnergy(const double& X, const double& Y, const double& Z) {
    return (X-Y+Z)/(2.*sqrt(X));
}
double DPhSp::ETemplateA(const double& X, const double& Y, const double& Z) {
    const double res = GetEnergy(X, Y, Z);
    if (res < 0) {
        cout << "ETemplateA: negative value " << res << endl;
        return -res;
    }
    return  res;
}
double DPhSp::ETemplateC(const double& X, const double& Y, const double& Z) {
    const double res = -GetEnergy(X, Y, Z);
    if (res < 0) {
        cout << "ETemplateC: negative value " << res << endl;
        return -res;
    }
    return res;
}
double DPhSp::ETemplateM(const double& X, const double& Y, const double& Z) {
    return  GetEnergy(X, Z, Y);
}
double DPhSp::Momentum(const double& E, const double& msq) {
    const double psq = E*E-msq;
    return (psq > 0 ? sqrt(psq) : 0);
}

double DPhSp::eA_AB(const double& mABsq, const double& mAsq,
                    const double& mBsq) {return ETemplateA(mABsq, mBsq, mAsq);}
double DPhSp::eB_AB(const double& mABsq, const double& mAsq,
                    const double& mBsq) {return ETemplateA(mABsq, mAsq, mBsq);}
double DPhSp::eC_AB(const double& mABsq, const double& mMsq,
                    const double& mCsq) {return ETemplateC(mABsq, mMsq, mCsq);}
double DPhSp::eM_AB(const double& mABsq, const double& mMsq,
                    const double& mCsq) {return ETemplateM(mABsq, mMsq, mCsq);}

double DPhSp::eA_AC(const double& mACsq, const double& mAsq,
                    const double& mCsq) {return ETemplateA(mACsq, mCsq, mAsq);}
double DPhSp::eB_AC(const double& mACsq, const double& mMsq,
                    const double& mBsq) {return ETemplateC(mACsq, mMsq, mBsq);}
double DPhSp::eC_AC(const double& mACsq, const double& mAsq,
                    const double& mCsq) {return ETemplateA(mACsq, mAsq, mCsq);}
double DPhSp::eM_AC(const double& mACsq, const double& mMsq,
                    const double& mBsq) {return ETemplateM(mACsq, mMsq, mBsq);}

double DPhSp::eA_BC(const double& mBCsq, const double& mMsq,
                    const double& mAsq) {return ETemplateC(mBCsq, mMsq, mAsq);}
double DPhSp::eB_BC(const double& mBCsq, const double& mBsq,
                    const double& mCsq) {return ETemplateA(mBCsq, mCsq, mBsq);}
double DPhSp::eC_BC(const double& mBCsq, const double& mBsq,
                    const double& mCsq) {return ETemplateA(mBCsq, mBsq, mCsq);}
double DPhSp::eM_BC(const double& mBCsq, const double& mMsq,
                    const double& mAsq) {return ETemplateM(mBCsq, mMsq, mAsq);}

double DPhSp::pA_AB(const double& mABsq, const double& mAsq,
                    const double& mBsq) {
    return Momentum(eA_AB(mABsq, mAsq, mBsq), mAsq);
}
double DPhSp::pB_AB(const double& mABsq, const double& mAsq,
                    const double& mBsq) {
    return Momentum(eB_AB(mABsq, mAsq, mBsq), mBsq);
}
double DPhSp::pC_AB(const double& mABsq, const double& mMsq,
                    const double& mCsq) {
    return Momentum(eC_AB(mABsq, mMsq, mCsq), mCsq);
}
double DPhSp::pM_AB(const double& mABsq, const double& mMsq,
                    const double& mCsq) {
    return Momentum(eM_AB(mABsq, mMsq, mCsq), mMsq);
}

double DPhSp::pA_AC(const double& mACsq, const double& mAsq,
                    const double& mCsq) {
    return Momentum(eA_AC(mACsq, mAsq, mCsq), mAsq);
}
double DPhSp::pB_AC(const double& mACsq, const double& mMsq,
                    const double& mBsq) {
    return Momentum(eB_AC(mACsq, mMsq, mBsq), mBsq);
}
double DPhSp::pC_AC(const double& mACsq, const double& mAsq,
                    const double& mCsq) {
    return Momentum(eC_AC(mACsq, mAsq, mCsq), mCsq);
}
double DPhSp::pM_AC(const double& mACsq, const double& mMsq,
                    const double& mBsq) {
    return Momentum(eM_AC(mACsq, mMsq, mBsq), mMsq);
}

double DPhSp::pA_BC(const double& mBCsq, const double& mMsq,
                    const double& mAsq) {
    return Momentum(eA_BC(mBCsq, mMsq, mAsq), mAsq);
}
double DPhSp::pB_BC(const double& mBCsq, const double& mBsq,
                    const double& mCsq) {
    return Momentum(eB_BC(mBCsq, mBsq, mCsq), mBsq);
}
double DPhSp::pC_BC(const double& mBCsq, const double& mBsq,
                    const double& mCsq) {
    return Momentum(eC_BC(mBCsq, mBsq, mCsq), mCsq);
}
double DPhSp::pM_BC(const double& mBCsq, const double& mMsq,
                    const double& mAsq) {
    return Momentum(eM_BC(mBCsq, mMsq, mAsq), mMsq);
}

double DPhSp::eAB_M(const double& mMsq, const double& mABsq,
                    const double& mCsq) {return ETemplateA(mMsq, mCsq, mABsq);}
double DPhSp::eAC_M(const double& mMsq, const double& mACsq,
                    const double& mBsq) {return ETemplateA(mMsq, mBsq, mACsq);}
double DPhSp::eBC_M(const double& mMsq, const double& mBCsq,
                    const double& mAsq) {return ETemplateA(mMsq, mAsq, mBCsq);}

double DPhSp::pAB_M(const double& mMsq, const double& mABsq,
                    const double& mCsq) {
    return Momentum(eAB_M(mMsq, mABsq, mCsq), mABsq);
}
double DPhSp::pAC_M(const double& mMsq, const double& mACsq,
                    const double& mBsq) {
    return Momentum(eAC_M(mMsq, mACsq, mBsq), mACsq);
}
double DPhSp::pBC_M(const double& mMsq, const double& mBCsq,
                    const double& mAsq) {
    return Momentum(eBC_M(mMsq, mBCsq, mAsq), mBCsq);
}

void DPhSp::GetLVs(const double& mABsq, const double& mACsq,
                   LVect* pM, LVect* pA, LVect* pB, LVect* pC) const {
    // Lorentz vectorz are calculater in M rest frame.
    // pA is directed to axis z
    // Decay plane is xz plane
    // pxB is chosen to be positive
    *pA = LVect((mABsq+mACsq-m_mBsq-m_mCsq)/(2.*m_mM), 0, 0, 0);
    *pB = LVect((m_mMsq+m_mBsq-mACsq)/(2.*m_mM), 0, 0, 0);
    *pC = LVect((m_mMsq+m_mCsq-mABsq)/(2.*m_mM), 0, 0, 0);
    *pM = LVect(m_mM, 0, 0, 0);

    pA->z(sqrt(std::fabs(pA->t()*pA->t() - m_mAsq)));
    pB->z((m_mBsq+m_mAsq+2*pA->t()*pB->t()-mABsq)/(2*pA->z()));
    pC->z((m_mCsq+m_mAsq+2*pA->t()*pC->t()-mACsq)/(2*pA->z()));

    pB->x(sqrt(std::fabs(pB->t()*pB->t()-pB->z()*pB->z()-m_mBsq)));
    pC->x(-pB->x());
}

