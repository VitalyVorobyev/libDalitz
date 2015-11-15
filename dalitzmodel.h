#ifndef DALITZMODEL_H
#define DALITZMODEL_H

#include "EvtResonance2.h"
#include "dalitzphasespace.h"
#include <vector>

class DalitzModel : public DalitzPhaseSpace{
public:
  DalitzModel(const double& mmo, const double& mcha, const double& mchb, const double& mchc);

  virtual EvtComplex Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3) = 0;

  EvtComplex Amp(const double& mAB, const double& mAC);
  double P(const double& mAB, const double& mAC);
  double Arg(const double& mAB, const double& mAC);

  void AddRes(EvtResonance2* res) {m_res_v.push_back(res); return;}
  const EvtResonance2* Res(const int resn) {return resn<m_res_v.size() ? m_res_v[resn] : m_res_v[0];}
  int ResNum(void) const {return m_res_v.size();}

  void SetGamma(const int resn, const double& a) {if(resn<m_res_v.size()) m_res_v[resn]->SetGamma(a); return;}
  void SetMass(const int resn, const double& a)  {if(resn<m_res_v.size()) m_res_v[resn]->SetMass(a);  return;}
  void SetAmp(const int resn, const double& a)   {if(resn<m_res_v.size()) m_res_v[resn]->SetAmp(a);   return;}
  void SetTheta(const int resn, const double& a) {if(resn<m_res_v.size()) m_res_v[resn]->SetTheta(a); return;}
  void SetMomenta(const int resn, const EvtVector4R& p4_p, const EvtVector4R& p4_d1,const EvtVector4R& p4_d2){if(resn<m_res_v.size()) m_res_v[resn]->SetMomenta(p4_p,p4_d1,p4_d2); return;}

  void SetABaxis(const std::string& str) {mABaxis = str; return;}
  void SetACaxis(const std::string& str) {mACaxis = str; return;}
  void SetBCaxis(const std::string& str) {mBCaxis = str; return;}

  std::string ABaxis(void) const { return mABaxis;}
  std::string ACaxis(void) const { return mACaxis;}
  std::string BCaxis(void) const { return mBCaxis;}
private:
  std::vector<EvtResonance2*> m_res_v;
  std::string mABaxis;
  std::string mACaxis;
  std::string mBCaxis;
};

#endif // DALITZMODEL_H
