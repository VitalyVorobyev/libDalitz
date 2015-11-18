#include "dalitzmcintegral.h"

DalitzMCIntegral::DalitzMCIntegral(const DalitzModel &_dm):
  RandomDalitzPoint(_dm),m_ncounts(1e6),m_int(0),m_int_err(0)
{
  dm = const_cast<DalitzModel*>(&_dm);
}

double DalitzMCIntegral::GetIntegral(const long& nc){
  if(nc) m_ncounts = nc;
  double res = 0;
  double ressq = 0;
  double mAB,mAC;
  std::cout << "DalitzMCIntegral: integrating..." << std::endl;
  for(int i=0; i<m_ncounts; i++){
    if(!(i%100000)) std::cout << i << " counts" << std::endl;
    GetPoint(mAB,mAC);
    res += dm->P(mAB,mAC);
    ressq += res*res;
  }
  m_int = res/m_ncounts;
  m_int_err = ressq/m_ncounts;
  std::cout << "Done! Int = " << m_int << " +- " << m_int_err << std::endl;
  return m_int;
}
