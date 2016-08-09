#ifndef ABSDALITZMODEL_H
#define ABSDALITZMODEL_H

#include <complex>
#include <string>
#include <vector>

#include "dalitzphasespace.h"
#include "dstrip.h"

///
/// \brief The MnPar class
///
class MnPar{
public:
  MnPar(const str& n, cdouble& iv, cdouble& ie, cdouble& rl, cdouble& ll):
    name(n),ival(iv),ierr(ie),rlim(rl),llim(ll)
  {}
  str name;
  double ival;
  double ierr;
  double rlim;
  double llim;
};

///
/// \brief The AbsDalitzModel class
///
class AbsDalitzModel : public DalitzPhaseSpace{
public:
  AbsDalitzModel(const DalitzPhaseSpace& phsp);
  AbsDalitzModel(cdouble& mM, cdouble& mA, cdouble& mB, cdouble& mC);

  /// Get DalitzModel complex amplitude
  compld Amp(cdouble& mABsq, cdouble& mACsq) const;
  /// Get modulo squared DalitzModel amplitude
  double P(cdouble& mABsq, cdouble& mACsq) const;
  /// Get complex phase (radians) of DalitzModel amplitude
  double Arg(cdouble& mABsq, cdouble& mACsq) const;

  /// Set caption for mAB^2 axis
  void SetABaxis(cstr& s) {mABaxis = s; return;}
  /// Set caption for mAC^2 axis
  void SetACaxis(cstr& s) {mACaxis = s; return;}
  /// Set caption for mBC^2 axis
  void SetBCaxis(cstr& s) {mBCaxis = s; return;}

  /// Get caption for mAB^2 axis
  str ABaxis(void) const { return mABaxis;}
  /// Get caption for mAC^2 axis
  str ACaxis(void) const { return mACaxis;}
  /// Get caption for mBC^2 axis
  str BCaxis(void) const { return mBCaxis;}

  /// Get name of the i'th DalitzPlotObject
  str ResName(const int resn) const {return m_res_names[resn];}
  /// Get name of the i'th DalitzPlotObject
  str AmpName(const int resn) const {return m_amp_names[resn];}

  /// To fill std::vector<compld> resv with amplitudes of each DalitzPlotObject of a DalitzModel. Returns size full amplitude
  compld GetAmplitudes(vectcd& resv, cdouble& mABsq, cdouble& mACsq) const;
  /// Dimention of the normalization matrix
  unsigned AmpNum(void) const { return m_ampl.size();}

  ///
  int OpenCachedIntegrals(const std::string& fname,const bool silent = true);
  ///
  /// \brief NormWithCache. Speed up the computation of normalization using the relation
  /// I = \sum_i |a_i|^2 I_i + 2Re(\sum_{i>j} a_i a_j^* I_{ij}), where
  /// I_i = \int |A_i|^2 dm_+^2 dm_-^2 and
  /// I_{ij} = \int A_i A_j^* dm_+^2 dm_-^2, where
  /// A_i is a complex amplitude of i^th resonance
  /// \return Value of normalization integral
  ///
  double NormWithCache(void) const;
  ///
  /// Vector of coef
  void GetCoefficients(vectcd& coefv) const;
  ///
  std::vector<DStrip*> GetResAreas(void) const {return m_res_areas;}
  ///
  const DStrip* GetResArea(const int i) const {return m_res_areas[i];}
  /// Get number of DalitzPlotObject in DalitzModel
  unsigned ResNum(void) const {return m_res_names.size();}
  /// Get vector of complex amplitudes for precalculation of the normalization integrals
  void GetAmpVals(vectcd& resv, cdouble& mABsq,cdouble& mACsq) const;
  /// Get
  void ShowAmpls(void) const;
  ///
  str GetAmpStr(void) const;

  /// Get vector of complex amplitudes for all resonances
  virtual void GetResVals(vectcd& resv, cdouble& mABsq,cdouble& mACsq) const = 0;
  /// Get
  virtual compld GetResVal(cdouble& mABsq,cdouble& mACsq, const int resnum) const = 0;
  /// Set parameters
  virtual void SetParams(const std::vector<double>& pars) = 0;
  /// Get normalization
  virtual double Norm(void) const = 0;
  /// Vector of parameter for minuit2
  std::vector<MnPar> MnPars(void) const;

protected:
  str mABaxis;
  str mACaxis;
  str mBCaxis;

  void SetResNames(    const vectstr& vnames) {m_res_names      = vnames;}
  void SetAmpNames(    const vectstr& vnames) {m_amp_names      = vnames;}
  void SetCoefficients(const vectcd&  coefv)  {m_ampl           = coefv;}
  void SetAmpSignature(const vecti&   vals)   {m_amp_signature  = vals;}
  void SetResAreas(const std::vector<DStrip*>& vec){m_res_areas = vec;}
  void SetResAreas(const vectd &ledge, const vectd &redge, const vecti& types);

private:
  /// Matrix of cached normalization integrals
  std::vector<vectcd> m_res_int;
  /// List of complex coefficients (Amplitudes) for normalization units
  vectcd m_ampl;
  /// Correspondence between Amplitudes and Resonances
  vecti m_amp_signature;
  /// List of resonance names
  vectstr m_res_names;
  /// List of normalization unit names
  vectstr m_amp_names;
  /// List of Dalitz Strips for normalization units
  std::vector<DStrip*> m_res_areas;
};

#endif // ABSDALITZMODEL_H
