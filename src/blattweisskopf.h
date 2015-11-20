/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtBlattWeisskopf.hh,v 1.2 2009-03-16 16:43:40 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

/// \brief Blatt-Weisskopf penetration form factor for a resonance R->AB.
/// Taken from CLEO preprint 00-23 (hep-ex/0011065)

#ifndef BLATTWEISSKOPF_H
#define BLATTWEISSKOPF_H

#include "formfactor.h"

class FFType{
public:
  static const int FFMeson     = 0;
  static const int FFResonance = 1;
};

class BlattWeisskopf : public FormFactor{
public:
  BlattWeisskopf(const int LL, const double& p0sq, const int type);
  BlattWeisskopf(const BlattWeisskopf&);
  ~BlattWeisskopf();

  double operator()(const double& psq) const;

  static double m_r_meson;     /// radius of a meson state
  static double m_r_resonance; /// radius of a resonance

private:
  int m_spin; /// angular momentum of resonance
  int m_type; /// FF type (meson or resonance)
  double m_F0;/// formula evaluated at p0

  double compute(const double& psq) const;
};

#endif // BLATTWEISSKOPF_H
