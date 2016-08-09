#include "dalitzphasespace.h"
#include <cmath>
#include <iostream>

using namespace std;

DalitzPhaseSpace::DalitzPhaseSpace(cdouble& mM, cdouble& mA, cdouble& mB, cdouble& mC):
  m_mM(mM),m_mA(mA),m_mB(mB),m_mC(mC),
  m_mMsq(m_mM*m_mM),m_mAsq(m_mA*m_mA),
  m_mBsq(m_mB*m_mB),m_mCsq(m_mC*m_mC),
  m_mABsq_min((m_mA+m_mB)*(m_mA+m_mB)),
  m_mABsq_max((m_mM-m_mC)*(m_mM-m_mC)),
  m_mACsq_min((m_mA+m_mC)*(m_mA+m_mC)),
  m_mACsq_max((m_mM-m_mB)*(m_mM-m_mB)),
  m_mBCsq_min((m_mB+m_mC)*(m_mB+m_mC)),
  m_mBCsq_max((m_mM-m_mA)*(m_mM-m_mA)),
  m_mass_sq_sum(m_mMsq+m_mAsq+m_mBsq+m_mCsq),
  m_area(0)
{
}

DalitzPhaseSpace::DalitzPhaseSpace(const DalitzPhaseSpace& phsp):
  DalitzPhaseSpace(phsp.mM(),phsp.mA(),phsp.mB(),phsp.mC())
{
}

double DPhSp::GetmBCsq(cdouble& mABsq,cdouble& mACsq) const {
  cdouble mBCsq = m_mass_sq_sum-mABsq-mACsq;
  if(mBCsq<0){
    cout << "DPhSp::GetmBCsq: " << m_mass_sq_sum << " " << mABsq << " " << mACsq << " " << mBCsq << endl;
    cout << "IsInPlot: " << IsInPlot(mABsq,mACsq) << endl;
  }
  return mBCsq;
}

/// Min value of AB squared invariant mass
double DPhSp::GetmABsqMin(cdouble& mA, cdouble& mB) {return pow(mA+mB,2);}

/// Max value of AB squared invariant mass
double DPhSp::GetmABsqMax(cdouble& mM, cdouble& mC) {return pow(mM-mC,2);}

// Particle energies in the resonance rest frame
double DPhSp::eA_AB(cdouble& mABsq) const {return eA_AB(mABsq,m_mAsq,m_mBsq);}
double DPhSp::eB_AB(cdouble& mABsq) const {return eB_AB(mABsq,m_mAsq,m_mBsq);}
double DPhSp::eC_AB(cdouble& mABsq) const {return eC_AB(mABsq,m_mMsq,m_mCsq);}
double DPhSp::eM_AB(cdouble& mABsq) const {return eM_AB(mABsq,m_mMsq,m_mCsq);}

double DPhSp::eA_AC(cdouble& mACsq) const {return eA_AC(mACsq,m_mAsq,m_mCsq);}
double DPhSp::eB_AC(cdouble& mACsq) const {return eB_AC(mACsq,m_mMsq,m_mBsq);}
double DPhSp::eC_AC(cdouble& mACsq) const {return eC_AC(mACsq,m_mAsq,m_mCsq);}
double DPhSp::eM_AC(cdouble& mACsq) const {return eM_AC(mACsq,m_mMsq,m_mBsq);}

double DPhSp::eA_BC(cdouble& mBCsq) const {return eA_BC(mBCsq,m_mMsq,m_mAsq);}
double DPhSp::eB_BC(cdouble& mBCsq) const {return eB_BC(mBCsq,m_mBsq,m_mCsq);}
double DPhSp::eC_BC(cdouble& mBCsq) const {return eC_BC(mBCsq,m_mBsq,m_mCsq);}
double DPhSp::eM_BC(cdouble& mBCsq) const {return eM_BC(mBCsq,m_mMsq,m_mAsq);}

// Particle momenta in the resonance rest frame
double DPhSp::pA_AB(cdouble& mABsq) const {return pA_AB(mABsq,m_mAsq,m_mBsq);}
double DPhSp::pB_AB(cdouble& mABsq) const {return pB_AB(mABsq,m_mAsq,m_mBsq);}
double DPhSp::pC_AB(cdouble& mABsq) const {return pC_AB(mABsq,m_mMsq,m_mCsq);}
double DPhSp::pM_AB(cdouble& mABsq) const {return pM_AB(mABsq,m_mMsq,m_mCsq);}

double DPhSp::pA_AC(cdouble& mACsq) const {return pA_AC(mACsq,m_mAsq,m_mCsq);}
double DPhSp::pB_AC(cdouble& mACsq) const {return pB_AC(mACsq,m_mMsq,m_mBsq);}
double DPhSp::pC_AC(cdouble& mACsq) const {return pC_AC(mACsq,m_mAsq,m_mCsq);}
double DPhSp::pM_AC(cdouble& mACsq) const {return pM_AC(mACsq,m_mMsq,m_mBsq);}

double DPhSp::pA_BC(cdouble& mBCsq) const {return pA_BC(mBCsq,m_mMsq,m_mAsq);}
double DPhSp::pB_BC(cdouble& mBCsq) const {return pB_BC(mBCsq,m_mBsq,m_mCsq);}
double DPhSp::pC_BC(cdouble& mBCsq) const {return pC_BC(mBCsq,m_mBsq,m_mCsq);}
double DPhSp::pM_BC(cdouble& mBCsq) const {return pM_BC(mBCsq,m_mMsq,m_mAsq);}

// Resonance energy in the M frame
double DPhSp::eAB_M(cdouble& mABsq) const {return eAB_M(m_mMsq,mABsq,m_mCsq);}
double DPhSp::eAC_M(cdouble& mACsq) const {return eAC_M(m_mMsq,mACsq,m_mBsq);}
double DPhSp::eBC_M(cdouble& mBCsq) const {return eBC_M(m_mMsq,mBCsq,m_mAsq);}

// Resonance momentum in the M frame
double DPhSp::pAB_M(cdouble& mABsq) const {return pAB_M(m_mMsq,mABsq,m_mCsq);}
double DPhSp::pAC_M(cdouble& mACsq) const {return pAC_M(m_mMsq,mACsq,m_mBsq);}
double DPhSp::pBC_M(cdouble& mBCsq) const {return pBC_M(m_mMsq,mBCsq,m_mAsq);}

// Helicity angles
double DPhSp::CosHelAB(cdouble& mABsq, cdouble& mACsq, double& pq) const {
  return CosHelAB(m_mMsq,m_mAsq,m_mBsq,m_mCsq,mABsq,mACsq,pq);
}
double DPhSp::CosHelAC(cdouble& mACsq, cdouble& mBCsq, double& pq) const {
  return CosHelAC(m_mMsq,m_mAsq,m_mBsq,m_mCsq,mACsq,mBCsq,pq);
}
double DPhSp::CosHelBC(cdouble& mBCsq, cdouble& mABsq, double& pq) const {
  return CosHelBC(m_mMsq,m_mAsq,m_mBsq,m_mCsq,mBCsq,mABsq,pq);
}

void DPhSp::RangeCalc(cdouble& e1,cdouble& e2,cdouble& m1sq,cdouble& m2sq,double& msqmin,double& msqmax){
  cdouble p1sq = e1*e1-m1sq;
  cdouble p1   = p1sq > 0 ? sqrt(p1sq) : 0;
  cdouble p2sq = e2*e2-m2sq;
  cdouble p2   = p2sq > 0 ? sqrt(p2sq) : 0;
  cdouble esum = e1+e2;
  cdouble pdif = p1-p2;
  cdouble psum = p1+p2;
  msqmax = (esum-pdif)*(esum+pdif);
  msqmin = (esum-psum)*(esum+psum);
}

int DPhSp::mABsqRange_AC(cdouble& mACsq, double& mABsqMin, double& mABsqMax) const{
  if(mACsq > m_mACsq_max || mACsq < m_mACsq_min) return -1;
  cdouble eA = eA_AC(mACsq);
  cdouble eB = eB_AC(mACsq);
//  cout << "eA: " << eA << ", eB: " << eB << ", mACsq: " << mACsq << endl;
  RangeCalc(eA,eB,m_mAsq,m_mBsq,mABsqMin,mABsqMax);
  return 0;
}

int DPhSp::mBCsqRange_AC(cdouble& mACsq, double& mBCsqMin, double& mBCsqMax) const{
  if(mACsq > m_mACsq_max || mACsq < m_mACsq_min) return -1;
  cdouble eC = eC_AC(mACsq);
  cdouble eB = eB_AC(mACsq);
  RangeCalc(eC,eB,m_mCsq,m_mBsq,mBCsqMin,mBCsqMax);
  return 0;
}

int DPhSp::mBCsqRange_AB(cdouble& mABsq, double& mBCsqMin, double& mBCsqMax) const{
  if(mABsq > m_mABsq_max || mABsq < m_mABsq_min) return -1;
  cdouble eB = eB_AB(mABsq);
  cdouble eC = eC_AB(mABsq);
  RangeCalc(eB,eC,m_mBsq,m_mCsq,mBCsqMin,mBCsqMax);
  return 0;
}

int DPhSp::mACsqRange_AB(cdouble& mABsq, double& mACsqMin, double& mACsqMax) const{
  if(mABsq > m_mABsq_max || mABsq < m_mABsq_min) return -1;
  cdouble eA = eA_AB(mABsq);
  cdouble eC = eC_AB(mABsq);
  RangeCalc(eA,eC,m_mAsq,m_mCsq,mACsqMin,mACsqMax);
  return 0;
}

int DPhSp::mACsqRange_BC(cdouble& mBCsq, double& mACsqMin, double& mACsqMax) const{
  if(mBCsq > m_mBCsq_max || mBCsq < m_mBCsq_min) return -1;
  cdouble eC = eC_BC(mBCsq);
  cdouble eA = eA_BC(mBCsq);
  RangeCalc(eC,eA,m_mCsq,m_mAsq,mACsqMin,mACsqMax);
  return 0;
}

int DPhSp::mABsqRange_BC(cdouble& mBCsq, double& mABsqMin, double& mABsqMax) const{
  if(mBCsq > m_mBCsq_max || mBCsq < m_mBCsq_min) return -1;
  cdouble eB = eB_BC(mBCsq);
  cdouble eA = eA_BC(mBCsq);
  RangeCalc(eB,eA,m_mBsq,m_mAsq,mABsqMin,mABsqMax);
  return 0;
}

bool DPhSp::IsInPlot(cdouble& mABsq,cdouble& mACsq) const {
  double mABsqmin,mABsqmax;
  if(mABsqRange_AC(mACsq,mABsqmin,mABsqmax)) return false;
  if(mABsq < mABsqmin || mABsq > mABsqmax)   return false;
//  cout << "IsInPlot: " << mABsq << " in (" << mABsqmin << "," << mABsqmax << ")" << endl;
  return true;
}

//double DPhSp::CosHelAB(cdouble& mABsq, cdouble& mBCsq, double& pq) const{
//  return CosHelAB(m_mMsq,m_mAsq,m_mBsq,m_mCsq,mABsq,mBCsq,pq);
//}

//double DPhSp::CosHelAC(cdouble& mACsq, cdouble& mBCsq, double& pq) const{
//  return CosHelAB(m_mMsq,m_mAsq,m_mBsq,m_mCsq,mABsq,mBCsq,pq);
//}

//double DPhSp::CosHelBC(cdouble& mBCsq, cdouble& mABsq, double& pq) const{
//  return CosHelAB(m_mMsq,m_mAsq,m_mBsq,m_mCsq,mABsq,mBCsq,pq);
//}

// * Static methods * //
double DPhSp::GetmBCsq(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq,cdouble& mACsq){
  return mMsq+mAsq+mBsq+mCsq-mABsq-mACsq;
}

double DPhSp::pRes(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq){
  cdouble res = (0.25*pow(mMsq-mAsq-mBsq,2)-mAsq*mBsq)/mMsq;
  return res>0 ? sqrt(res) : 0;
}

double DPhSp::q(cdouble& mRsq, cdouble& mA, cdouble& mB){
  cdouble var = (mRsq-pow(mA+mB,2))*(mRsq-pow(mA-mB,2));
  return var>0 ? sqrt(var/(4.*mRsq)) : 0;
}

double DPhSp::p(cdouble& mRsq, cdouble& mM, cdouble& mC){
  cdouble var = (pow(mM-mC,2)-mRsq)*(pow(mM+mC,2)-mRsq);
  return var>0 ? sqrt(var/(4.*mRsq)) : 0;
}

double DPhSp::ysq(cdouble& mMsq, cdouble& mABsq, cdouble& mCsq){
// ??? to be checked
    return 0.5*(mMsq + mABsq - mCsq)/sqrt(mABsq*mMsq) - 1.;
}

double DPhSp::CosHelAB(cdouble& mMo,cdouble& mA,cdouble& mB,cdouble& mC,cdouble& mABsq){
  cdouble mABmin = GetmABsqMin(mA,mB);
  cdouble mABmax = GetmABsqMax(mMo,mC);
  return (mABmax+mABmin-2.*mABsq)/(mABmax-mABmin);
}

//double DPhSp::CosHelAB(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq, cdouble& mBCsq, double& pq){
//  cdouble enB   = eB_AB(mABsq,mAsq,mBsq);
//  cdouble moBsq = pow(enB,2) - mBsq;
//  if(moBsq<0){ pq = 0; return 0;}
//  cdouble moB   = sqrt(moBsq);
//  cdouble enC   = eC_AB(mABsq,mMsq,mCsq);
//  cdouble moCsq = pow(enC,2) - msq;
//  if(moCsq<0){ pq = 0; return 0;}
//  cdouble moC   = sqrt(moCsq);
//  pq = moB*moC;
//  return (mBsq + mCsq + 2.*enB*enC - mBCsq)/(2.*moB*moC);
//}

double DPhSp::CosHelAB(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq, cdouble& mACsq, double& pq){
  cdouble eA   = eA_AB(mABsq,mAsq,mBsq);
  cdouble pAsq = pow(eA,2) - mAsq;
  if(pAsq<0){ pq = 0; return 0;}
  cdouble pA   = sqrt(pAsq);
  cdouble eC   = eC_AB(mABsq,mMsq,mCsq);
  cdouble pCsq = pow(eC,2) - mCsq;
  if(pCsq<0){ pq = 0; return 0;}
  cdouble pC   = sqrt(pCsq);
  pq = pA*pC;
  return (mAsq + mCsq + 2.*eA*eC - mACsq)/(2.*pA*pC);
}

double DPhSp::CosHelAC(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mACsq, cdouble& mBCsq, double& pq){
  cdouble eB   = eB_AC(mACsq,mMsq,mBsq);
  cdouble pBsq = pow(eB,2) - mBsq;
  if(pBsq<0){ pq = 0; return 0;}
  cdouble pB   = sqrt(pBsq);
  cdouble eC   = eC_AC(mACsq,mAsq,mCsq);
  cdouble pCsq = pow(eC,2) - mCsq;
  if(pCsq<0){ pq = 0; return 0;}
  cdouble pC   = sqrt(pCsq);
  pq = pB*pC;
  return (mBsq + mCsq + 2.*eB*eC - mBCsq)/(2.*pB*pC);
}

double DPhSp::CosHelBC(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mBCsq, cdouble& mABsq, double& pq){
  cdouble eA   = eA_BC(mBCsq,mMsq,mAsq);
  cdouble pAsq = pow(eA,2) - mAsq;
  if(pAsq<0){ pq = 0; return 0;}
  cdouble pA   = sqrt(pAsq);
  cdouble eB   = eB_BC(mBCsq,mBsq,mCsq);
  cdouble pBsq = pow(eB,2) - mBsq;
  if(pBsq<0){ pq = 0; return 0;}
  cdouble pB   = sqrt(pBsq);
  pq = pA*pB;
  return (mAsq + mBsq + 2.*eA*eB - mABsq)/(2.*pA*pB);
}

//double DPhSp::ETemplateA(cdouble& X,cdouble& Y,cdouble& Z){return (X-Y+Z)/sqrt(2.*X);}
//double DPhSp::ETemplateC(cdouble& X,cdouble& Y,cdouble& Z){return (Y-X-Z)/sqrt(2.*X);}
//double DPhSp::ETemplateM(cdouble& X,cdouble& Y,cdouble& Z){return (Y+X-Z)/sqrt(2.*X);}
double DPhSp::GetEnergy( cdouble& X,cdouble& Y,cdouble& Z){return (X-Y+Z)/(2.*sqrt(X));}
double DPhSp::ETemplateA(cdouble& X,cdouble& Y,cdouble& Z){
  cdouble res = GetEnergy(X,Y,Z);
  if(res<0){
    cout << "ETemplateA: negative value " << res << endl;
    return -res;
  }
  return  res;
}
double DPhSp::ETemplateC(cdouble& X,cdouble& Y,cdouble& Z){
  cdouble res = -GetEnergy(X,Y,Z);
  if(res<0){
    cout << "ETemplateC: negative value " << res << endl;
    return -res;
  }
  return res;
}
double DPhSp::ETemplateM(cdouble& X,cdouble& Y,cdouble& Z){return  GetEnergy(X,Z,Y);}
double DPhSp::Momentum(cdouble& E, cdouble& msq){
  cdouble psq = E*E-msq;
  return (psq > 0 ? sqrt(psq) : 0);
}

double DPhSp::eA_AB(cdouble& mABsq,cdouble& mAsq,cdouble& mBsq) {return ETemplateA(mABsq,mBsq,mAsq);}
double DPhSp::eB_AB(cdouble& mABsq,cdouble& mAsq,cdouble& mBsq) {return ETemplateA(mABsq,mAsq,mBsq);}
double DPhSp::eC_AB(cdouble& mABsq,cdouble& mMsq,cdouble& mCsq) {return ETemplateC(mABsq,mMsq,mCsq);}
double DPhSp::eM_AB(cdouble& mABsq,cdouble& mMsq,cdouble& mCsq) {return ETemplateM(mABsq,mMsq,mCsq);}

double DPhSp::eA_AC(cdouble& mACsq,cdouble& mAsq,cdouble& mCsq) {return ETemplateA(mACsq,mCsq,mAsq);}
double DPhSp::eB_AC(cdouble& mACsq,cdouble& mMsq,cdouble& mBsq) {return ETemplateC(mACsq,mMsq,mBsq);}
double DPhSp::eC_AC(cdouble& mACsq,cdouble& mAsq,cdouble& mCsq) {return ETemplateA(mACsq,mAsq,mCsq);}
double DPhSp::eM_AC(cdouble& mACsq,cdouble& mMsq,cdouble& mBsq) {return ETemplateM(mACsq,mMsq,mBsq);}

double DPhSp::eA_BC(cdouble& mBCsq,cdouble& mMsq,cdouble& mAsq) {return ETemplateC(mBCsq,mMsq,mAsq);}
double DPhSp::eB_BC(cdouble& mBCsq,cdouble& mBsq,cdouble& mCsq) {return ETemplateA(mBCsq,mCsq,mBsq);}
double DPhSp::eC_BC(cdouble& mBCsq,cdouble& mBsq,cdouble& mCsq) {return ETemplateA(mBCsq,mBsq,mCsq);}
double DPhSp::eM_BC(cdouble& mBCsq,cdouble& mMsq,cdouble& mAsq) {return ETemplateM(mBCsq,mMsq,mAsq);}

double DPhSp::pA_AB(cdouble& mABsq,cdouble& mAsq,cdouble& mBsq) {return Momentum(eA_AB(mABsq,mAsq,mBsq),mAsq);}
double DPhSp::pB_AB(cdouble& mABsq,cdouble& mAsq,cdouble& mBsq) {return Momentum(eB_AB(mABsq,mAsq,mBsq),mBsq);}
double DPhSp::pC_AB(cdouble& mABsq,cdouble& mMsq,cdouble& mCsq) {return Momentum(eC_AB(mABsq,mMsq,mCsq),mCsq);}
double DPhSp::pM_AB(cdouble& mABsq,cdouble& mMsq,cdouble& mCsq) {return Momentum(eM_AB(mABsq,mMsq,mCsq),mMsq);}

double DPhSp::pA_AC(cdouble& mACsq,cdouble& mAsq,cdouble& mCsq) {return Momentum(eA_AC(mACsq,mAsq,mCsq),mAsq);}
double DPhSp::pB_AC(cdouble& mACsq,cdouble& mMsq,cdouble& mBsq) {return Momentum(eB_AC(mACsq,mMsq,mBsq),mBsq);}
double DPhSp::pC_AC(cdouble& mACsq,cdouble& mAsq,cdouble& mCsq) {return Momentum(eC_AC(mACsq,mAsq,mCsq),mCsq);}
double DPhSp::pM_AC(cdouble& mACsq,cdouble& mMsq,cdouble& mBsq) {return Momentum(eM_AC(mACsq,mMsq,mBsq),mMsq);}

double DPhSp::pA_BC(cdouble& mBCsq,cdouble& mMsq,cdouble& mAsq) {return Momentum(eA_BC(mBCsq,mMsq,mAsq),mAsq);}
double DPhSp::pB_BC(cdouble& mBCsq,cdouble& mBsq,cdouble& mCsq) {return Momentum(eB_BC(mBCsq,mBsq,mCsq),mBsq);}
double DPhSp::pC_BC(cdouble& mBCsq,cdouble& mBsq,cdouble& mCsq) {return Momentum(eC_BC(mBCsq,mBsq,mCsq),mCsq);}
double DPhSp::pM_BC(cdouble& mBCsq,cdouble& mMsq,cdouble& mAsq) {return Momentum(eM_BC(mBCsq,mMsq,mAsq),mMsq);}

double DPhSp::eAB_M(cdouble& mMsq, cdouble& mABsq, cdouble& mCsq) {return ETemplateA(mMsq,mCsq,mABsq);}
double DPhSp::eAC_M(cdouble& mMsq, cdouble& mACsq, cdouble& mBsq) {return ETemplateA(mMsq,mBsq,mACsq);}
double DPhSp::eBC_M(cdouble& mMsq, cdouble& mBCsq, cdouble& mAsq) {return ETemplateA(mMsq,mAsq,mBCsq);}

double DPhSp::pAB_M(cdouble& mMsq, cdouble& mABsq, cdouble& mCsq) {return Momentum(eAB_M(mMsq,mABsq,mCsq),mABsq);}
double DPhSp::pAC_M(cdouble& mMsq, cdouble& mACsq, cdouble& mBsq) {return Momentum(eAC_M(mMsq,mACsq,mBsq),mACsq);}
double DPhSp::pBC_M(cdouble& mMsq, cdouble& mBCsq, cdouble& mAsq) {return Momentum(eBC_M(mMsq,mBCsq,mAsq),mBCsq);}

void DPhSp::GetLVs(cdouble& mABsq,cdouble& mACsq, EvtVector4R& pM, EvtVector4R& pA,EvtVector4R& pB,EvtVector4R& pC) const{
  // Lorentz vectorz are calculater in M rest frame.
  // pA is directed to axis z
  // Decay plane is xz plane
  // pxB is chosen to be positive
  pA = EvtVector4R((mABsq+mACsq-m_mBsq-m_mCsq)/(2.*m_mM),0,0,0);
  pB = EvtVector4R((m_mMsq+m_mBsq-mACsq)/(2.*m_mM),0,0,0);
  pC = EvtVector4R((m_mMsq+m_mCsq-mABsq)/(2.*m_mM),0,0,0);
  pM = EvtVector4R(m_mM,0,0,0);

  pA.pz(sqrt(fabs(pA.e()*pA.e() - m_mAsq)));
  pB.pz((m_mBsq+m_mAsq+2*pA.e()*pB.e()-mABsq)/(2*pA.pz()));
  pC.pz((m_mCsq+m_mAsq+2*pA.e()*pC.e()-mACsq)/(2*pA.pz()));

  pB.px(sqrt(fabs(pB.e()*pB.e()-pB.pz()*pB.pz()-m_mBsq)));
  pC.px(-pB.px());
  return;
}

