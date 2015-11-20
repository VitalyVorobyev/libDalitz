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

class DalitzResonance : public DalitzPlotObject{
public:
  DalitzResonance(const std::string& name,const int PropType,const int WidthType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const EvtComplex& camp);               /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const int WidthType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const double& amp, const double& phi); /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const EvtComplex& camp);               /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const double& amp, const double& phi); /// Constuctor for RBW and GS resonances
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& alpha,const EvtComplex& camp); /// Constuctor for NR
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& alpha,const double& amp,const double& phi); /// Constuctor for NR
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const EvtComplex& camp); /// Constuctor for Bugg f0(500)
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& amp,const double& phi); /// Constuctor for Bugg f0(500)
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& beta1,const double& beta2, const EvtComplex& amp); /// Constructor for virtual D*(2010)
  DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& beta1,const double& beta2, const double& amp,const double& phi); /// Constructor for virtual D*(2010)

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
