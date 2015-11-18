#ifndef DALITZPLOTOBJECT_H
#define DALITZPLOTOBJECT_H

#include <string>
#include "EvtComplex.h"

class DalitzPlotObject{
public:
  DalitzPlotObject(const std::string& name, const EvtComplex& amp = EvtComplex(0,0));
  DalitzPlotObject(const std::string& name, const double& a, const double& phi);
  virtual EvtComplex evaluate() = 0;

  void SetName(const std::string& name) {m_name = name; return;}
  std::string Name(void) const {return m_name;}

  void SetCAmp(const EvtComplex& amp);
  EvtComplex CAmp(void) const {return m_amp;}

  void SetAmp(const double& a);
  double Amp(void) const {return m_a;}

  void SetPhase(const double& phi);
  double Phase(void) const {return m_phi;}
private:
  std::string m_name;
  EvtComplex m_amp;
  double m_a;
  double m_phi;
};

#endif // DALITZPLOTOBJECT_H
