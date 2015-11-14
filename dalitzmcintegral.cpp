#include "dalitzmcintegral.h"

DalitzMCIntegral::DalitzMCIntegral(const DalitzModel &_dm):
  RandomDalitzPoint(_dm),m_ncounts(1e6)
{
  dm = const_cast<DalitzModel*>(&_dm);
}

double DalitzMCIntegral::GetIntegral(const long& nc){
  if(nc) m_ncounts = nc;
  double res = 0;
  double mAB,mAC;
  for(int i=0; i<m_ncounts; i++){
    GetPoint(mAB,mAC);
    res += dm->P(mAB,mAC);
  }
  return res/m_ncounts;
}
