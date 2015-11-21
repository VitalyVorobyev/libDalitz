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
  EvtComplex amp(0.,0.);
  for(int i=0; i<(int)m_res_v.size(); i++) amp += GetResAmp(m_res_v[i],mAB,mAC);
  return amp;
}

EvtComplex DalitzModel::GetResAmp(const DalitzPlotObject* res,const double& mAB, const double& mAC){
  switch(res->Path()){
  case ResPath::AB: return res->evaluate(mAC,mBCsq(mAC,mAB));
  case ResPath::AC: return res->evaluate(mAB,mBCsq(mAC,mAB));
  case ResPath::BC: return res->evaluate(mAC,mAB);
  }
  return EvtComplex(0.,0.);
}

EvtComplex DalitzModel::GetAmplitudes(std::vector<EvtComplex>& vec,const double& mAB, const double& mAC){
  vec.clear();
  EvtComplex amp(0.,0.);
  for(int i=0; i<(int)m_res_v.size(); i++){
    vec.push_back(GetResAmp(m_res_v[i],mAB,mAC));
    amp += vec[i];
  }
  return amp;
}
