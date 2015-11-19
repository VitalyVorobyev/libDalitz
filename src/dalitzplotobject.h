#ifndef DALITZPLOTOBJECT_H
#define DALITZPLOTOBJECT_H

#include <string>
#include "EvtComplex.h"

/// Abstract class for any object can appears on a Dalitz diagram.
/// The class contains name, complex amplitude of the object and virtual
/// method evaluate()

class DalitzPlotObject{
public:
  DalitzPlotObject(const std::string& name, const EvtComplex& amp = EvtComplex(0,0));
  DalitzPlotObject(const std::string& name, const double& a, const double& phi);
  virtual EvtComplex evaluate(const double& mACsq,const double& mBCsq) = 0;

  /// Modificators
  void SetName(const std::string& name) {m_name = name; return;}
  void SetCAmp(const EvtComplex& amp);
  void SetAmp(const double& a);
  void SetPhase(const double& phi);

  /// Interface
  std::string Name(void) const {return m_name;}
  EvtComplex CAmp(void) const {return m_amp;}
  double Amp(void) const {return m_a;}
  double Phase(void) const {return m_phi;}

private:
  std::string m_name;
  EvtComplex m_amp;
  double m_a;
  double m_phi;
};

#endif // DALITZPLOTOBJECT_H
