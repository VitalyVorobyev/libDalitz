#ifndef BTODK_DTOKSPIPI_MODEL_H
#define BTODK_DTOKSPIPI_MODEL_H

#include "kspipimodel.h"

///
/// \brief The BtoDK_DtoKspipi_Model class. The class describes CP-violating amplitude
/// of the B- -> D0 K-, D0 -> Ks0 pi+ pi- decay.
///

class BtoDK_DtoKspipi_Model : public KspipiModel{
public:
  BtoDK_DtoKspipi_Model(const double& gamma = 70,const double& delb = 120, const double& rb = 0.1);

  compld Amp(const double& mAB, const double& mAC) const;
  /// Interface for Minuit2 FCN
  void SetParams(const std::vector<double>& par);

  void SetFlv(const int flv);
  void SetGD(const double& gamma, const double& delta, const double& rB);
  void SetXY(const double& xp, const double& yp,const double& xm, const double& ym);

  int Flv(void)      const;
  double rBp(void)   const;
  double rBm(void)   const;
  double xp(void)    const;
  double xm(void)    const;
  double yp(void)    const;
  double ym(void)    const;
  double gamma(void) const;
  double delb(void)  const;

private:
  int SetParamsGD(const std::vector<double>& par);
  int SetParamsXY(const std::vector<double>& par);

  compld m_zp;
  compld m_zm;
  int m_flv;
};

#endif // BTODK_DTOKSPIPI_MODEL_H

