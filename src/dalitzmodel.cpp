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
  EvtVector4R pM;
  EvtVector4R pA;
  EvtVector4R pB;
  EvtVector4R pC;
  GetLVs(mAB,mAC,pM,pA,pB,pC);
  return Amp(pM,pA,pB,pC);
}
