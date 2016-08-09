#ifndef DALITZRESONANCE_H
#define DALITZRESONANCE_H

#include "dalitzplotobject.h"
#include "resdecayangulardistribution.h"
#include "formfactor.h"
#include "abspropagator.h"
#include "dalitzphasespace.h"
#include "consts.h"

#include <string>

///
/// \brief The ResPath class. Class is keeping flags specifying resonance path (R->AB, R->AC or R->BC)
///
class ResPath{
public:
  static const int AB = 0;
  static const int AC = 1;
  static const int BC = 2;
};

///
/// \brief The DalitzResonance class. Class describes complex amplitude of a three-body decay through an intermediate resonance (M -> RC, R -> AB)
///
class DalitzResonance : public DalitzPlotObject{
public:
  /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const int WidthType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const compld& camp);
  /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const int WidthType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const double& amp, const double& phi);
  /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const compld& camp);
  /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const double& amp, const double& phi);
  /// Constuctor for NR and Flatte
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& alpha,const compld& camp);
  /// Constuctor for NR and Flatte
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& alpha,const double& amp,const double& phi);
  /// Constuctor for Bugg f0(500) and virtual D*0(2010) by Belle
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const compld& camp);
  /// Constuctor for Bugg f0(500) and virtual D*0(2010) by Belle
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& amp,const double& phi);
  /// Constructor for virtual D*(2010) by LHCb and rho-omega interference
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& beta1,const double& beta2, const compld& amp);
  /// Constructor for virtual D*(2010) by LHCb and rho-omega interference
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& beta1,const double& beta2, const double& amp,const double& phi);

  compld evaluate(const double& mACsq, const double& mBCsq) const;

  int Path(void) const {return m_path;}

  /// Set if the resonance amplitude was changed
  void SetAmpUpd(const bool x) {m_amp_upd = x;}
  /// Set if the resonance parameters were changed
  void SetParUpd(const bool x) {m_par_upd = x;}
  /// Check if the resonance amplitude was changed
  bool IsAmpUpd(void) const {return m_amp_upd;}
  /// Check if the resonance parameters were changed
  bool IsParUpd(void) const {return m_par_upd;}

private:
  double SetFFAngAmp(const DalitzPhaseSpace* phsp, const double &mres, const bool vdst = false);/// Initializes form factors and angular distribution

  ResDecayAngularDistribution* m_ang_amp;
  FormFactor* m_mff;
  FormFactor* m_rff;
  AbsPropagator* m_prop;

  int m_ptype;
  int m_wtype;
  int m_spin;
  int m_path;

  bool m_amp_upd;
  bool m_par_upd;
};

#endif // DALITZRESONANCE_H
