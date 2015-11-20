#include "dalitzmodel.h"

DalitzModel::DalitzModel(const double& mmo,const double& mcha,const double& mchb,const double& mchc) :
DalitzPhaseSpace(mmo,mcha,mchb,mchc)
{
}

double DalitzModel::P(const double& mAB,const double& mAC){
  return abs2(Amp(mAB,mAC));
}

double DalitzModel::Arg(const double& mAB,const double& mAC){
  return arg(Amp(mAB,mAC));
}

EvtComplex DalitzModel::Amp(const double& mAB, const double& mAC){
  EvtComplex res(0.,0.);
  for(int i=0; i<(int)m_res_v.size(); i++){
    switch(m_res_v[i]->Path()){
    case ResPath::AB:
      res += m_res_v[i]->evaluate(mAC,mBCsq(mAC,mAB));
//      res += m_res_v[i]->evaluate(mBC(mAC,mAB),mAC);
      break;
    case ResPath::AC:
      res += m_res_v[i]->evaluate(mAB,mBCsq(mAC,mAB));
//      res += m_res_v[i]->evaluate(mBC(mAC,mAB),mAB);
      break;
    case ResPath::BC:
//      res += m_res_v[i]->evaluate(mAB,mAC);
      res += m_res_v[i]->evaluate(mAC,mAB);
      break;
    default:
      break;
    }
  }
  return res;// + EvtComplex(-2.537,0.923);
}
