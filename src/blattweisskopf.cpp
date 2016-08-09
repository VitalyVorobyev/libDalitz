#include "blattweisskopf.h"
#include <math.h>
#include <iostream>

double BlattWeisskopf::m_r_meson     = 5.0;// GeV-1
double BlattWeisskopf::m_r_resonance = 1.5;// GeV-1

BlattWeisskopf::BlattWeisskopf(const int spin,const double& p0, const int type):
  FormFactor(type == FFType::FFMeson ? m_r_meson : m_r_resonance,p0),m_spin(spin),
  m_type(type),m_F0(compute(p0))
{}

BlattWeisskopf::BlattWeisskopf(const BlattWeisskopf& other):
   FormFactor(other.r(),other.p0()),m_spin(other.m_spin),m_F0(other.m_F0)
{}

BlattWeisskopf::~BlattWeisskopf()
{}

double BlattWeisskopf::operator()(const double& p) const {
  return compute(p)/m_F0;
}

// Blatt-Weisskopf form factors
// see e.g. hep-ex/0011065
// Dalitz Analysis of the Decay D0->K-pi+pi0 (CLEO)
//
// p   - momentum of either daugher in the meson rest frame,
//       the mass of the meson is used
// pAB - momentum of either daughter in the candidate rest frame
//       the mass of the candidate is used
// R   - meson radial parameter
//
// In the CLEO paper R=5 GeV-1 for D0, R=1.5 for intermediate resonances

double BlattWeisskopf::compute(const double& p) const {
  if(!m_spin) return 1.;
  double denom(1.0);
  const double z = p*r();
  const double zSq = z*z;

  switch(m_spin){
  case 0:
    denom = 1.0;
    break;
  case 1:
    denom = 1.0 + zSq;
    break;
  case 2:
    denom = zSq*(zSq + 3.0) + 9.0;
    break;
  case 3:
    denom = zSq*(zSq*(zSq + 6.0) + 45.0) + 225.0;
    break;
  case 4:
    denom = zSq*(zSq*(zSq*(zSq + 10.0) + 135.0) + 1575.0) + 11025.0;
    break;
  case 5:
    denom = zSq*(zSq*(zSq*(zSq*(zSq + 15.0) + 315.0) + 6300.0) + 99225.0) + 893025.0;
    break;
  default:
    std::cout << "BlattWeisskopf::compute: wrong spin " << m_spin << std::endl;
    break;
  }
  return sqrt(1.0/denom);
}
