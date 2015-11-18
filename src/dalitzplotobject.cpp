#include "dalitzplotobject.h"

DalitzPlotObject::DalitzPlotObject(const std::string& name, const EvtComplex& amp):
  m_name(name)
{
  SetCAmp(amp);
}

DalitzPlotObject::DalitzPlotObject(const std::string& name, const double& a, const double& phi):
  DalitzPlotObject(name,a*EvtComplex(cos(phi),sin(phi)))
{
}

void DalitzPlotObject::SetPhase(const double& phi) {
  m_phi = phi;
  m_amp = m_a*EvtComplex(cos(m_phi),sin(m_phi));
  return;
}

void DalitzPlotObject::SetAmp(const double& a) {
  m_a = a;
  m_amp = m_a*EvtComplex(cos(m_phi),sin(m_phi));
  return;
}

void DalitzPlotObject::SetCAmp(const EvtComplex& amp){
  m_amp = amp;
  m_a = abs(amp); m_phi = arg(amp);
  return;
}
