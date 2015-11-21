#include "resdecayangulardistribution.h"
#include "dalitzphasespace.h"
#include <math.h>

ResDecayAngularDistribution::ResDecayAngularDistribution(const int spin, const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres):
  m_spin(spin),m_mmo_sq(mmo*mmo),m_mca_sq(mca*mca),m_mcb_sq(mcb*mcb),m_mcc_sq(mcc*mcc),m_mre_sq(mres*mres)
{
       if(m_spin == 1) make_precalc_one();
  else if(m_spin == 2) make_precalc_all();
}

double ResDecayAngularDistribution::operator()(const double& mACsq,const double& mBCsq) const{
  double pdf = 1;
  switch(m_spin){
  case 0:
    pdf = 1.;
    break;
  case 1:
    pdf = mACsq-mBCsq+m_pc_var1/m_mre_sq;
    break;
  case 2:
    const double mABsq = DalitzPhaseSpace::mBCsq(m_mmo_sq,m_mca_sq,m_mcb_sq,m_mcc_sq,mACsq,mBCsq);
    pdf = pow(mACsq-mBCsq+m_pc_var1/m_mre_sq,2)-
          (1.0/3.0)*(mABsq-m_pc_var2+m_pc_var3/m_mre_sq)*
                    (mABsq-m_pc_var4+m_pc_var5/m_mre_sq);
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
