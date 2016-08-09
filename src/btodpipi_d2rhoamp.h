#ifndef BTODPIPI_D2RHOAMP_H
#define BTODPIPI_D2RHOAMP_H

#include <vector>
#include <complex>
#include <string>

#include "absdalitzmodel.h"

class BtoDpipi_D2rhoAmp : public AbsDalitzModel{
public:
  BtoDpipi_D2rhoAmp();
  /// Vector of Res_amp (rho, omega and rho' combined)
  void GetResVals(vectcd& resv,cdouble& mABsq,cdouble& mACsq) const;
  /// Complex amplitude of a resonance
  compld GetResVal(cdouble& mABsq,cdouble& mACsq, const int resnum) const;
  /// Set parameters
  void SetParams(const vectd& pars);
  /// Get normalization
  double Norm(void) const;

private:
  void init(void);

  compld AmpABspin2(cdouble& mABsq,cdouble& mR,cdouble& gR) const;
  compld AmpBCspin1(cdouble& mBCsq,cdouble& mR,cdouble& gR) const;

  compld DdstAmp(cdouble& mABsq,cdouble& helAB,cdouble& pqAB) const;
  compld rhoAmp( cdouble& mBCsq,cdouble& helBC,cdouble& pqBC) const;

  double hwrho(cdouble& s, double& hwd) const;

  double ang1(cdouble& hel,cdouble& pq) const;
  double ang2(cdouble& hel,cdouble& pq) const;

  double BWUnit1(cdouble& p, cdouble& scale) const;
  double BWUnit2(cdouble& p, cdouble& scale) const;
  double BlattWeisskopf(cdouble& p0,cdouble& p, cint spin, cdouble& scale) const;

  double m_ffB;
  double m_ffR;

  double m_Ddst_mass;
  double m_Ddst_width;
  double m_rho_mass;
  double m_rho_width;

  compld m_Ddst_amp;
  compld m_rho_amp;
};

#endif // BTODPIPI_D2RHOAMP_H
