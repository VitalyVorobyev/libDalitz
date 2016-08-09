#include "resdecayangulardistribution.h"
#include "dalitzphasespace.h"
#include <math.h>

bool ResDecayAngularDistribution::m_use_mRsq = true;

ResDecayAngularDistribution::ResDecayAngularDistribution(const int spin, const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres):
  m_spin(spin),m_mmo(mmo),m_mca(mca),m_mcb(mcb),m_mcc(mcc),
  m_mmo_sq(mmo*mmo),m_mca_sq(mca*mca),m_mcb_sq(mcb*mcb),m_mcc_sq(mcc*mcc),m_mre_sq(mres*mres)
{
       if(m_spin == 1) make_precalc_one();
  else if(m_spin == 2) make_precalc_all();
}

double ResDecayAngularDistribution::operator()(const double& mACsq,const double& mBCsq) const{
  if(!m_spin) return 1.;
  double pdf = 1;
  const double mABsq = DalitzPhaseSpace::GetmBCsq(m_mmo_sq,m_mca_sq,m_mcb_sq,m_mcc_sq,mACsq,mBCsq);
//  double ysq     = m_spin>2 ? DalitzPhaseSpace::ysq(m_mmo_sq,mABsq,m_mcc_sq) : 1;
  double cos_hel = m_spin>2 ? DalitzPhaseSpace::CosHelAB(m_mmo,m_mca,m_mcb,m_mcc,mABsq) : 1;
  double p       = m_spin>2 ? DalitzPhaseSpace::p(m_mre_sq,m_mmo,m_mcc) : 0;
  double q       = m_spin>2 ? DalitzPhaseSpace::q(m_mre_sq,m_mca,m_mcb) : 0;
  double mRsq = m_mre_sq;
  switch(m_spin){
  case 1:
    if(!m_use_mRsq) mRsq = mABsq;
    pdf = mACsq-mBCsq+m_pc_var1/mRsq;
    break;
  case 2:
    if(!m_use_mRsq) mRsq = mABsq;
    pdf = pow(mACsq-mBCsq+m_pc_var1/mRsq,2)-
          (1.0/3.0)*(mABsq-m_pc_var2+m_pc_var3/mRsq)*
                    (mABsq-m_pc_var4+m_pc_var5/mRsq);
    break;
  case 3:// spin 3 resonance Zemach helicity factor
//    pdf = sqrt(1.+ysq)*(1.+0.4*ysq)*(pow(cos_hel,3)-0.6*cos_hel)*pow(p*q,3);
    pdf = -8.0*(5.0*pow(cos_hel,3) - 3.0*cos_hel)/5.0*pow(p*q,3);
    break;
  case 4:// spin 4 resonance Zemach helicity factor
//    pdf = ((8.*ysq*ysq+40.*ysq)/35.+1)*(pow(cos_hel,4)-(30.*cos_hel*cos_hel+3.)/35.)*pow(p*q,4);
    pdf = 16.0*(35.0*pow(cos_hel,4) - 30.0*pow(cos_hel,2) + 3.0)/35.0*pow(p*q,4);
    break;
  case 5:// spin 5 resonance Zemach helicity factor
    pdf = -32.0*(63.0*pow(cos_hel,5) - 70.0*pow(cos_hel,3) + 15.0*cos_hel)/63.0*pow(p*q,5);
  }
// ** Laura++ code ** //
  // Calculate the spin factors
  // These are calculated as follows
  // -2^j * (q*p)^j * cj * Pj(cosHel)
  // where Pj(coshHel) is the jth order Legendre polynomial and
  // cj = j! / (2j-1)!!

//  Double_t spinTerm(1.0);
//  if (resSpin_ == 1) {
//          // Calculate vector resonance Zemach helicity factor
//          spinTerm = -2.0*q_*p_*cosHel;
//  } else if (resSpin_ == 2) {
//          // Calculate tensor resonance Zemach helicity factor
//          Double_t pProd = q_*p_;
//          spinTerm = 4.0*(pProd*pProd)*(3.0*cosHel*cosHel - 1.0)/3.0;
//  } else if (resSpin_ == 3) {
//          // Calculate spin 3 resonance Zemach helicity factor
//          Double_t pProd = q_*p_;
//          spinTerm = -8.0*(pProd*pProd*pProd)*(5.0*cosHel*cosHel*cosHel - 3.0*cosHel)/5.0;
//  } else if (resSpin_ == 4) {
//          // Calculate spin 4 resonance Zemach helicity factor
//          Double_t pProd = q_*p_;
//          spinTerm = 16.0*(pProd*pProd*pProd*pProd)*(35.0*cosHel*cosHel*cosHel*cosHel - 30.0*cosHel*cosHel + 3.0)/35.0;
//  } else if (resSpin_ == 5) {
//          // Calculate spin 5 resonance Zemach helicity factor
//          Double_t pProd = q_*p_;
//          spinTerm = -32.0*(pProd*pProd*pProd*pProd*pProd)*(63.0*cosHel*cosHel*cosHel*cosHel*cosHel - 70.0*cosHel*cosHel*cosHel + 15.0*cosHel)/63.0;
//  }
// ** ** //

//  std::cout << "ResDecayAngAmp: " << pdf << std::endl;
  return pdf;
}

void ResDecayAngularDistribution::make_precalc_one(void){
  m_pc_var1 = (m_mmo_sq - m_mcc_sq)*(m_mcb_sq - m_mca_sq);
  return;
}

void ResDecayAngularDistribution::make_precalc_all(void){
  make_precalc_one();
  m_pc_var2 = 2.*(m_mmo_sq+m_mcc_sq);
  m_pc_var3 = pow(m_mmo_sq-m_mcc_sq,2);
  m_pc_var4 = 2.*(m_mca_sq+m_mcb_sq);
  m_pc_var5 = pow(m_mca_sq-m_mcb_sq,2);
}
