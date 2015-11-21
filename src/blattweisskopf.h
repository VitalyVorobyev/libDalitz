/// \brief Blatt-Weisskopf penetration form factor for a resonance R->AB.
/// Taken from CLEO preprint 00-23 (hep-ex/0011065)
/// See original paper J. Blatt and V. Weisskopf, "Theoretical Nuclear Physics"
/// (Wiley, New York, 1952).

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
  BlattWeisskopf(const int spin, const double& p0, const int type);
  BlattWeisskopf(const BlattWeisskopf&);
  ~BlattWeisskopf();

  double operator()(const double& p) const;

  static double m_r_meson;     /// radius of a meson state
  static double m_r_resonance; /// radius of a resonance

private:
  int m_spin; /// angular momentum of resonance
  int m_type; /// FF type (meson or resonance)
  double m_F0;/// formula evaluated at p0

  double compute(const double& psq) const;
};

#endif // BLATTWEISSKOPF_H
