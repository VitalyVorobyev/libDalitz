/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtBlattWeisskopf.hh,v 1.2 2009-03-16 16:43:40 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Blatt-Weisskopf penetration form factor for a resonance R->AB.
// Taken from CLEO preprint 00-23 (hep-ex/0011065)

#ifndef BLATTWEISSKOPF_H
#define BLATTWEISSKOPF_H

#include "formfactor.h"

class BlattWeisskopf : public FormFactor{
public:
  BlattWeisskopf(const int LL, const double& R, const double& _p0);
  BlattWeisskopf(const BlattWeisskopf&);
  ~BlattWeisskopf();

  double operator()(const double& p) const;

private:
  int _LL;// angular momentum of daughters
  double _F0; // formula evaluated at _p0
  double compute(const double& p) const;
};

#endif // BLATTWEISSKOPF_H
