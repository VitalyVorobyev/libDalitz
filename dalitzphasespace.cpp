#include "dalitzphasespace.h"

DalitzPhaseSpace::DalitzPhaseSpace(const double& mM, const double& mA, const double& mB, const double& mC):
m_mMo(mM),m_mChA(mA),m_mChB(mB),m_mChC(mC)
{
  m_mMo_sq  = m_mMo*m_mMo;
  m_mChA_sq = m_mChA*m_mChA;
  m_mChB_sq = m_mChB*m_mChB;
  m_mChC_sq = m_mChC*m_mChC;
  m_mAB_min = (m_mChA+m_mChB)*(m_mChA+m_mChB);
  m_mAC_min = (m_mChA+m_mChC)*(m_mChA+m_mChC);
  m_mBC_min = (m_mChB+m_mChC)*(m_mChB+m_mChC);
  m_mAB_max = (m_mMo-m_mChC)*(m_mMo-m_mChC);
  m_mAC_max = (m_mMo-m_mChB)*(m_mMo-m_mChB);
  m_mBC_max = (m_mMo-m_mChA)*(m_mMo-m_mChA);
}

DalitzPhaseSpace::DalitzPhaseSpace(const DalitzPhaseSpace& phsp):
  DalitzPhaseSpace(phsp.mM(),phsp.mA(),phsp.mB(),phsp.mC())
{
}

double DalitzPhaseSpace::mBC(const double& mAB,const double& mAC){
  return m_mMo_sq+m_mChA_sq+m_mChB_sq+m_mChC_sq-mAB-mAC;
}

// a particle energy in a resonance rest frame
//double DalitzPhaseSpace::eA_AB(const double& mAB) const{return 0.5*(mAB-m_mChB_sq+m_mChA_sq)/sqrt(mAB);}
//double DalitzPhaseSpace::eB_AB(const double& mAB) const{return 0.5*(mAB-m_mChA_sq+m_mChB_sq)/sqrt(mAB);}
//inline double DalitzPhaseSpace::eC_AB(const double& mAB){return 0.5*(m_mMo_sq- mAB-m_mChC_sq)/sqrt(mAB);}
inline double DalitzPhaseSpace::eA_AC(const double& mAC) const {return 0.5*(mAC-m_mChC_sq+m_mChA_sq)/sqrt(mAC);}
inline double DalitzPhaseSpace::eB_AC(const double& mAC) const {return 0.5*(m_mMo_sq- mAC-m_mChB_sq)/sqrt(mAC);}
//inline double DalitzPhaseSpace::eC_AC(const double& mAC){return 0.5*(mAC-m_mChA_sq+m_mChC_sq)/sqrt(mAC);}

int DalitzPhaseSpace::mAB_range(const double& mAC, double& mABmin, double& mABmax) const{
  if(mAC > m_mAC_max || mAC < m_mAC_min) return -1;
  const double eA  = eA_AC(mAC);
  const double eB  = eB_AC(mAC);
  const double pA  = sqrt(eA*eA-m_mChA_sq);
  const double pB  = sqrt(eB*eB-m_mChB_sq);
  const double esum = eA+eB;
  const double pdif = pA-pB;
  const double psum = pA+pB;
  mABmax = (esum-pdif)*(esum+pdif);
  mABmin = (esum-psum)*(esum+psum);
  return 0;
}

bool DalitzPhaseSpace::IsInPlot(const double& mAB,const double& mAC){
  double mABmin,mABmax;
  if(mAB_range(mAC,mABmin,mABmax)) return false;
  if(mAB >= mABmin && mAB <= mABmax) return true;
  else return false;
}

void DalitzPhaseSpace::GetLVs(const double& mAB,const double& mAC, EvtVector4R& pM, EvtVector4R& pA,EvtVector4R& pB,EvtVector4R& pC){
  // Lorentz vectorz are calculater in M rest frame.
  // pA is directed to axis z
  // Decay plane is xz plane
  // pxB is chosen to be positive
  const double mM_sq = m_mMo*m_mMo;
  const double mA_sq = m_mChA*m_mChA;
  const double mB_sq = m_mChB*m_mChB;
  const double mC_sq = m_mChC*m_mChC;

  double eA, pxA = 0 ,pyA = 0, pzA;
  double eB, pxB, pyB = 0, pzB;
  double eC, pxC, pyC = 0, pzC;
  double eM = m_mMo, pxM = 0,pyM = 0,pzM = 0;
  double mAB_test, mAC_test;

  eA  = (mAB+mAC-mB_sq-mC_sq)/(2.*m_mMo);
  eB  = (mM_sq+mB_sq-mAC)/(2.*m_mMo);
  eC  = (mM_sq+mC_sq-mAB)/(2.*m_mMo);

  pzA = sqrt(eA*eA - mA_sq);
  pzB = (mB_sq+mA_sq+2*eA*eB-mAB)/(2*pzA);
  pzC = (mC_sq+mA_sq+2*eA*eC-mAC)/(2*pzA);

  const double pxBsq = eB*eB - pzB*pzB - mB_sq;
  pxB = pxBsq > 0 ? sqrt(pxBsq) : 0;
  pxC = -pxB;

  mAC_test = (eA+eC)*(eA+eC)-(pxA+pxC)*(pxA+pxC)-(pyA+pyC)*(pyA+pyC)-(pzA+pzC)*(pzA+pzC);
  mAB_test = (eA+eB)*(eA+eB)-(pxA+pxB)*(pxA+pxB)-(pyA+pyB)*(pyA+pyB)-(pzA+pzB)*(pzA+pzB);

  if(fabs(mAC-mAC_test)>0.0001 || fabs(mAB-mAB_test)>0.0001 || std::isnan(pxB)){
    const double mB_test = sqrt(eB*eB-pxB*pxB-pyB*pyB-pzB*pzB);
    const double mC_test = sqrt(eC*eC-pxC*pxC-pyC*pyC-pzC*pzC);
    const double mA_test = sqrt(eA*eA-pxA*pxA-pyA*pyA-pzA*pzA);
    cout << "Wrong (mAC,mAB): (" << mAC << "," << mAB << ") -> (" << mAC_test << "," << mAB_test << "):" << endl;
    cout << "Masses squared: " << mM_sq << ", " << mA_sq << ", " << mB_sq << ", " << mC_sq << endl;
    cout << " B: (" << eB << "," << pxB << "," << pyB << "," << pzB << ") -> " << mB_test << endl;
    cout << " C: (" << eC << "," << pxC << "," << pyC << "," << pzC << ") -> " << mC_test << endl;
    cout << " A: (" << eA << "," << pxA << "," << pyA << "," << pzA << ") -> " << mA_test << endl;
  }
  pM = EvtVector4R(eM,pxM,pyM,pzM);
  pA = EvtVector4R(eA,pxA,pyA,pzA);
  pB = EvtVector4R(eB,pxB,pyB,pzB);
  pC = EvtVector4R(eC,pxC,pyC,pzC);
  return;
}
