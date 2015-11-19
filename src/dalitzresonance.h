#ifndef DALITZRESONANCE_H
#define DALITZRESONANCE_H

#include "dalitzplotobject.h"
#include "resdecayangulardistribution.h"
//#include "absvarwidth.h"
#include "formfactor.h"
#include "abspropagator.h"

#include <string>

/// \brief Class describes complex amplitude of a three-body decay through an intermediate resonance (M -> RC, R -> AB)
/// Class describes complex amplitude of a three-body decay through an intermediate resonance (M -> RC, R -> AB)

class ResPropType{
public:
  static const int NR;
  static const int RBW;
  static const int GS;
  static const int RhoOmega;
  static const int Bugg;
  static const int VDst;
};

class VarWType{
public:
  static const int Const;
  static const int BW;
  static const int GS;
  static const int Flatte;
  static const int Bugg;
};

class DalitzResonance : public DalitzPlotObject{
public:
  DalitzResonance(const std::string& name,const int PropType,const int WidthType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const EvtComplex& camp);               /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const int WidthType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const double& amp, const double& phi); /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const EvtComplex& camp);               /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const double& amp, const double& phi); /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const double& alpha,const EvtComplex& camp); /// Constuctor for NR
  DalitzResonance(const std::string& name,const int PropType,const double& alpha,const double& amp,const double& phi); /// Constuctor for NR
  DalitzResonance(const std::string& name,const int PropType,const EvtComplex& camp); /// Constuctor for Bugg f0(500)
  DalitzResonance(const std::string& name,const int PropType,const double& amp,const double& phi); /// Constuctor for Bugg f0(500)
//  DalitzResonance(const std::string& name,const int PropType,const double& amp, const double& phi); /// Constuctor for NR and Bugg
//  DalitzResonance(const std::string& name,const int PropType,const int WidthType, const double& mres, const EvtComplex& camp);               /// Constuctor for Flatte
//  DalitzResonance(const std::string& name,const int PropType,const int WidthType, const double& mres, const double& amp, const double& phi); /// Constuctor for Flatte

  EvtComplex evaluate(const double& mACsq,const double& mBCsq);

private:
  ResDecayAngularDistribution* m_ang_amp;
  FormFactor* m_mff;
  FormFactor* m_rff;
  AbsPropagator* m_prop;

  int m_ptype;
  int m_wtype;

//  double m_mot_sq;
//  double m_mca_sq;
//  double m_mcb_sq;
//  double m_mcc_sq;
};

#endif // DALITZRESONANCE_H