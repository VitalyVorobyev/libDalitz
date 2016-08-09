#include "dalitzmodel.h"
#include "dalitzmcintegral.h"

#include <iostream>

using namespace std;

DalitzModel::DalitzModel(cdouble& mmo,cdouble& mcha,cdouble& mchb,cdouble& mchc) :
AbsDalitzModel(mmo,mcha,mchb,mchc),use_subset(false)
{
}

bool DalitzModel::IsVetoed(cdouble& mAB, cdouble& mAC) const{
  for(unsigned i=0; i<m_veto_v.size(); i++){if((*m_veto_v[i])(mAB,mAC)) return true;}
  return false;
}

compld DalitzModel::Amp(cdouble& mAB, cdouble& mAC) const{
  compld amp(0.,0.);
  if(IsVetoed(mAB,mAC)) return amp;
  if(use_subset){// If some subset is chosen
    for(int i=0; i<(int)m_rlist.size(); i++) amp += GetResAmp(m_res_v[m_rlist[i]],mAB,mAC);
  } else{
    for(int i=0; i<(int)m_res_v.size(); i++) amp += GetResAmp(m_res_v[i],mAB,mAC);
  }
  return amp;
}

void DalitzModel::GetResVals(vectcd& resv, cdouble& mABsq,cdouble& mACsq) const{
  resv.clear();
  if(IsVetoed(mABsq,mACsq)){
    resv.resize(m_res_v.size(),0);
    return;
  }
  if(use_subset){// If some subset is chosen
    for(int i=0; i<(int)m_rlist.size(); i++) resv.push_back(GetResAmp(m_res_v[m_rlist[i]],mABsq,mACsq));
  } else{
    for(int i=0; i<(int)m_res_v.size(); i++) resv.push_back(GetResAmp(m_res_v[i],mABsq,mACsq));
  }
}

//vector<MnPar> DalitzModel::MnPars(void){
//  vector<MnPar> parv;
//  for(auto res : m_res_v){
//    parv.push_back(MnPar(res->Amp(),res->Amp()*0.1,res->Amp()*0.1,res->Amp()*10.));
//    parv.push_back(MnPar(res->Phase(),0.1,-2.*M_PI,2.*M_PI));
//  }
//  return parv;
//}

compld DalitzModel::ResAmp(const unsigned n, cdouble& mAB, cdouble& mAC) const{
  if(n<0 || n>=m_res_v.size()) return 0;
  return GetResAmp(m_res_v[n],mAB,mAC);
}

compld DalitzModel::GetResAmp(const DalitzPlotObject* res,cdouble& mAB, cdouble& mAC) const {
  switch(res->Path()){
  case ResPath::AB: return res->evaluate(mAC,GetmBCsq(mAC,mAB));
  case ResPath::AC: return res->evaluate(mAB,GetmBCsq(mAC,mAB));
  case ResPath::BC: return res->evaluate(mAC,mAB);
//  case ResPath::BC: return res->evaluate(mAB,mAC);
  }
  return compld(0.,0.);
}

compld DalitzModel::GetAmplitudes(vectcd& vec,cdouble& mAB, cdouble& mAC) const {
  vec.clear();
  compld amp(0.,0.);
  for(int i=0; i<(int)m_res_v.size(); i++){
    vec.push_back(GetResAmp(m_res_v[i],mAB,mAC));
    amp += vec[i];
  }
  return amp;
}

bool DalitzModel::SetRVec(const vector<int>& rlist){
  m_rlist.clear(); use_subset = false;
  cout << "Resonances subset: ";
  for(int i=0; i<(int)rlist.size(); i++){
    if(rlist[i]<(int)m_res_v.size()){
      cout << rlist[i] << " " << m_res_v[rlist[i]]->Name() << ", ";
      const int index = rlist[i];
      m_rlist.push_back(index);
    } else{
      cout << endl << "Wrong resonance index " << rlist[i] << " (should be less than " << m_res_v.size() << ")" << endl;
    }
  }
  cout << endl;
  use_subset = true;
  return true;
}

int DalitzModel::Norm(double& val, double& err, const long nc){
  DalitzMCIntegral mcint(this);
  return mcint.CalcIntegral(val,err,nc);
}

