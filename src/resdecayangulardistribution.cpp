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
  double ysq = 0;
  double cos_hel = 1;
  double p = 0;
  double q = 0;
  const double mABsq = DalitzPhaseSpace::mBCsq(m_mmo_sq,m_mca_sq,m_mcb_sq,m_mcc_sq,mACsq,mBCsq);
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
  case 3:
    ysq     = DalitzPhaseSpace::ysq(m_mmo_sq,mABsq,m_mcc_sq);
    cos_hel = DalitzPhaseSpace::CosHelAB(m_mmo,m_mca,m_mcb,m_mcc,mABsq);
    p       = DalitzPhaseSpace::p(m_mre_sq,m_mmo,m_mcc);
    q       = DalitzPhaseSpace::q(m_mre_sq,m_mca,m_mcb);
    pdf = sqrt(1.+ysq)*(1.+0.4*ysq)*(pow(cos_hel,3)-0.6*cos_hel)*pow(p*q,3);
    break;
  case 4:
    ysq     = DalitzPhaseSpace::ysq(m_mmo_sq,mABsq,m_mcc_sq);
    cos_hel = DalitzPhaseSpace::CosHelAB(m_mmo,m_mca,m_mcb,m_mcc,mABsq);
    p       = DalitzPhaseSpace::p(m_mre_sq,m_mmo,m_mcc);
    q       = DalitzPhaseSpace::q(m_mre_sq,m_mca,m_mcb);
    pdf = ((8.*ysq*ysq+40.*ysq)/35.+1)*(pow(cos_hel,4)-(30.*cos_hel*cos_hel+3.)/35.)*pow(p*q,4);
    break;
  }
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
