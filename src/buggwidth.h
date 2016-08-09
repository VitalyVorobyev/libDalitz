#ifndef BUGGWIDTH_H
#define BUGGWIDTH_H

#include "absvarwidth.h"

class BuggWidth : public AbsVarWidth{
public:
  BuggWidth(void);
//  ~BuggWidth() {}
  double mrGamma1(const double& s);
  void GetWidths(const double& s, double& G1, double& GTot);
  double sA(void) const {return m_sA;}
  double mrsq(void) const {return m_mrsq;}
  double g1sq(void) const {return m_g1sq_pc;}
  double z(void) const {return m_z_pc;}

  double operator()(const double& s = 0, const double& p = 0) const;

private:
  double m_mr,m_mrsq;
  double m_sA;
  double m_b1;
  double m_b2;
  double m_A;
  double m_g4pi;
  double m_alpha;

  double m_g1sq_pc;
  double m_rho1_ps;
  double m_z_pc;

  inline double rho_4pi(const double& s);
  inline double rho(const double& m, const double& s);

  double j1(const double& s);
  inline double z(const double& s);
  inline double g1sq(const double& s);

  double mrGamma2(const double& s);
  double mrGamma3(const double& s);
  double mrGamma4(const double& s);
};

#endif // BUGGWIDTH_H
