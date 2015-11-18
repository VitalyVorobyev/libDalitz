#ifndef FLATTEWIDTH_H
#define FLATTEWIDTH_H

#include "absvarwidth.h"

class FlatteWidth : public AbsVarWidth{
public:
  FlatteWidth(const double& m);

  double operator()(const double& s, const double& p) const;

private:
  double m_g1;
  double m_g2;
  double m_pi_sq;
  double m_pi0_sq;
  double m_K_sq;
  double m_K0_sq;

  double rho_pipi(const double& s) const;
  double rho_KK(const double& s) const;
  double phsp_factor(const double& msq, const double&s) const;
};

#endif // FLATTEWIDTH_H
