#ifndef DALITZRESONANCE_H
#define DALITZRESONANCE_H

#include "dalitzplotobject.h"
#include "resdecayangulardistribution.h"
#include "formfactor.h"
#include "abspropagator.h"
#include "dalitzphasespace.h"

#include <string>

/// \brief Class is keeping flags specifying resonance path (R->AB, R->AC or R->BC)

class ResPath{
public:
  static const int AB = 0;
  static const int AC = 1;
  static const int BC = 2;
};

/// \brief Class describes complex amplitude of a three-body decay through an intermediate resonance (M -> RC, R -> AB)

class DalitzResonance : public DalitzPlotObject{
public:
  DalitzResonance(const std::string& name,const int PropType,const int WidthType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const EvtComplex& camp);               /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const int WidthType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const double& amp, const double& phi); /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const EvtComplex& camp);               /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const double& amp, const double& phi); /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& alpha,const EvtComplex& camp);              /// Constuctor for NR and Flatte
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& alpha,const double& amp,const double& phi); /// Constuctor for NR and Flatte
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const EvtComplex& camp);              /// Constuctor for Bugg f0(500)
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& amp,const double& phi); /// Constuctor for Bugg f0(500)
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& beta1,const double& beta2, const EvtComplex& amp);               /// Constructor for virtual D*(2010) and rho-omega interference
  DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& beta1,const double& beta2, const double& amp,const double& phi); /// Constructor for virtual D*(2010) and rho-omega interference

  EvtComplex evaluate(const double& mACsq,const double& mBCsq);

  int Path(void) const {return m_path;}

private:
  double SetFFAngAmp(const DalitzPhaseSpace* phsp, const double &mres);/// Initializes form factors and angular distribution

  ResDecayAngularDistribution* m_ang_amp;
  FormFactor* m_mff;
  FormFactor* m_rff;
  AbsPropagator* m_prop;

  int m_ptype;
  int m_wtype;
  int m_spin;
  int m_path;
};

#endif // DALITZRESONANCE_H
