/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzplotobject.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/dalitzplotobject.h"

#include <cmath>

typedef std::complex<double> compld;

using std::sin;
using std::cos;

DalitzPlotObject::DalitzPlotObject(const std::string& name,
                                   const compld& amp) : m_name(name) {
    SetCAmp(amp);
}

DalitzPlotObject::DalitzPlotObject(const std::string& name,
                                   const double& a, const double& phi) :
    DalitzPlotObject(name, a*compld(cos(phi), sin(phi))) {}

void DalitzPlotObject::SetPhase(const double& phi) {
    m_phi = phi;
    m_amp = m_a*compld(cos(m_phi), sin(m_phi));
}

void DalitzPlotObject::SetAmp(const double& a) {
    m_a = a;
    m_amp = m_a*compld(cos(m_phi), sin(m_phi));
}

void DalitzPlotObject::SetCAmp(const compld& amp) {
    m_amp = amp;
    m_a = std::abs(amp); m_phi = std::arg(amp);
}
