#ifndef DALITZMODEL_H
#define DALITZMODEL_H

#include "dalitzphasespace.h"
#include "dalitzplotobject.h"
#include "dalitzresonance.h"

#include <vector>
#include <string>

/// \brief Class for full description of a three-body decay model
/// A DalitzModel class object contains std::vector of objects inheritanced
/// from DalitzPlotObject class. Method EvtComplex Amp(const double& mAB, const double& mAC)
/// returns coherent sum of complex amplitudes of DalitzPlotObject's.

class DalitzModel : public DalitzPhaseSpace{
public:
  DalitzModel(const double& mmo, const double& mcha, const double& mchb, const double& mchc);

  EvtComplex Amp(const double& mAB, const double& mAC); /// Get DalitzModel complex amplitude
  double P(const double& mAB, const double& mAC);       /// Get modulo squared DalitzModel amplitude
  double Arg(const double& mAB, const double& mAC);     /// Get complex phase (radians) of DalitzModel amplitude

  EvtComplex GetAmplitudes(std::vector<EvtComplex>& vec, const double& mAB, const double& mAC); /// To fill std::vector<EvtComplex> vec with amplitudes of each DalitzPlotObject of a DalitzModel. Returns size nomber of DalitzPlotObjects in the DalitzModel
  void AddRes(DalitzPlotObject* res) {m_res_v.push_back(res); return;} /// To add DalitzResonance (or DalitzPlotObject) to a DalitzModel
  const DalitzPlotObject* Res(const int resn) {return resn<(int)m_res_v.size() ? m_res_v[resn] : m_res_v[0];} /// Get pointer to i'th DalitzPlotObject of DalitzModel
  int ResNum(void) const {return m_res_v.size();} /// Get number of DalitzPlotObject in DalitzModel

//  void SetGamma(const int resn, const double& a) {if(resn<(int)m_res_v.size()) m_res_v[resn]->SetGamma(a); return;}
//  void SetMass(const int resn, const double& a)  {if(resn<(int)m_res_v.size()) m_res_v[resn]->SetMass(a);  return;}
  void SetAmp(const int resn, const double& a)   {if(resn<(int)m_res_v.size()) m_res_v[resn]->SetAmp(a);   return;} /// Set amplitude modulus of i'th DalitzPlotObject of DalitzModel
  void SetPhase(const int resn, const double& a) {if(resn<(int)m_res_v.size()) m_res_v[resn]->SetPhase(a); return;} /// Set amplitude phase of i'th DalitzPlotObject of DalitzModel
//  void SetMomenta(const int resn, const EvtVector4R& p4_p, const EvtVector4R& p4_d1,const EvtVector4R& p4_d2){
//    if(resn<(int)m_res_v.size()) m_res_v[resn]->SetMomenta(p4_p,p4_d1,p4_d2);
//    return;
//  }

  void SetABaxis(const std::string& str) {mABaxis = str; return;} /// Set caption for mAB^2 axis
  void SetACaxis(const std::string& str) {mACaxis = str; return;} /// Set caption for mAC^2 axis
  void SetBCaxis(const std::string& str) {mBCaxis = str; return;} /// Set caption for mBC^2 axis

  std::string ABaxis(void) const { return mABaxis;} /// Get caption for mAB^2 axis
  std::string ACaxis(void) const { return mACaxis;} /// Get caption for mAC^2 axis
  std::string BCaxis(void) const { return mBCaxis;} /// Get caption for mBC^2 axis

private:
  EvtComplex GetResAmp(const DalitzPlotObject *res, const double& mAB, const double& mAC); /// Chooses the correct combination of invariant masses squared and calculates resonance amplitude
  std::vector<DalitzPlotObject*> m_res_v;
  std::string mABaxis;
  std::string mACaxis;
  std::string mBCaxis;
};

#endif // DALITZMODEL_H
