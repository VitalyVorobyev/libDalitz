#ifndef DALITZMODEL_H
#define DALITZMODEL_H

#include "absdalitzmodel.h"

#include "dalitzplotobject.h"
#include "dalitzresonance.h"
#include "dalitzveto.h"
#include "consts.h"

#include <vector>
#include <string>

///
/// \brief Class for full description of a three-body decay model
/// A DalitzModel class object contains std::vector of objects inheritanced
/// from DalitzPlotObject class. Method Amp(const double& mAB, const double& mAC)
/// returns coherent sum of complex amplitudes of DalitzPlotObject's.
///
class DalitzModel : public AbsDalitzModel{
public:
  DalitzModel(const double& mmo, const double& mcha, const double& mchb, const double& mchc);

  /// Get DalitzModel complex amplitude
  compld Amp(const double& mAB, const double& mAC) const;

  /// Get resonance complex amplitude
  compld ResAmp(const unsigned n, const double& mAB, const double& mAC) const;

  /// Check if the Dalitz plot point is vetoed
  bool IsVetoed(const double& mAB, const double& mAC) const;

  /// To fill std::vector<compld> vec with amplitudes of each DalitzPlotObject of a DalitzModel. Returns size full amplitude
  compld GetAmplitudes(std::vector<compld>& vec, const double& mAB, const double& mAC) const ;

  /// To add DalitzResonance (or DalitzPlotObject) to a DalitzModel
  void AddRes(DalitzPlotObject* res) {m_res_v.push_back(res); return;}

  /// Get pointer to i'th DalitzPlotObject of DalitzModel
  const DalitzPlotObject* Res(const int resn) {return resn<(int)m_res_v.size() ? m_res_v[resn] : m_res_v[0];}

  ///
  double Norm(void) const {return 1;}

  ///
  void SetParams(const std::vector<double>& pars) {}

  /// Get vector of complex amplitudes for precalculation of the normalization integrals
  void GetResVals(vectcd& resv, cdouble& mABsq,cdouble& mACsq) const;
  /// Get
  compld GetResVal(cdouble& mABsq,cdouble& mACsq, const int resnum) const {return ResAmp(resnum,mABsq,mACsq);}

  /// To add DalitzVeto to a DalitzModel
  void AddVeto(DalitzVeto* veto) {m_veto_v.push_back(veto); return;}

//  void SetGamma(const int resn, const double& a) {if(resn<(int)m_res_v.size()) m_res_v[resn]->SetGamma(a); return;}
//  void SetMass(const int resn, const double& a)  {if(resn<(int)m_res_v.size()) m_res_v[resn]->SetMass(a);  return;}
  /// Set amplitude modulus of the i'th DalitzPlotObject of DalitzModel
  void SetAmp(const int resn, const double& a)   {if(resn<(int)m_res_v.size()) m_res_v[resn]->SetAmp(a);   return;}

  /// Set amplitude phase of the i'th DalitzPlotObject of DalitzModel
  void SetPhase(const int resn, const double& a) {if(resn<(int)m_res_v.size()) m_res_v[resn]->SetPhase(a); return;}
//  void SetMomenta(const int resn, const EvtVector4R& p4_p, const EvtVector4R& p4_d1,const EvtVector4R& p4_d2){
//    if(resn<(int)m_res_v.size()) m_res_v[resn]->SetMomenta(p4_p,p4_d1,p4_d2);
//    return;
//  }

  /// Set the list of resonances to be involved in the amplitude calculation
  bool SetRVec(const std::vector<int>& rlist);
  /// Calculate Monte Carlo normalization integral
  int Norm(double& val, double& err, const long nc = 0);

  int UpdateNormMap(void);
  double NormElement(const int i,const int j) const;

//  std::vector<MnPar> MnPars(void);

private:
  /// Chooses the correct combination of invariant masses squared and calculates resonance amplitude
  compld GetResAmp(const DalitzPlotObject *res, const double& mAB, const double& mAC) const ;
  std::vector<DalitzPlotObject*> m_res_v;
  std::vector<DalitzVeto*> m_veto_v;
  int use_subset;
  std::vector<int> m_rlist;
  std::vector<DStrip*> m_res_areas;

  std::vector<std::vector<double>> m_norm_map;
};

#endif // DALITZMODEL_H
