/** Copyright 2017 Vitaly Vorobyev
 ** @file btodk_dtokspipi_model.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_BTODK_DTOKSPIPI_MODEL_H_
#define SRC_BTODK_DTOKSPIPI_MODEL_H_

#include <complex>
#include <vector>

#include "./kspipimodel.h"

/**
 * @brief The BtoDK_DtoKspipi_Model class. The class describes CP-violating amplitude
 * of the B- -> D0 K-, D0 -> Ks0 pi+ pi- decay.
 */
class BtoDK_DtoKspipi_Model : public KspipiModel {
 public:
    BtoDK_DtoKspipi_Model(const double& gamma = 70, const double& delb = 120,
                          const double& rb = 0.1);
    std::complex<double> Amp(const double& mAB, const double& mAC) const;
    /// Interface for Minuit2 FCN
    void SetParams(const std::vector<double>& par);

    void SetFlv(const int flv);
    void SetGD(const double& gamma, const double& delta, const double& rB);
    void SetXY(const double& xp, const double& yp,
               const double& xm, const double& ym);

    int Flv(void) const;
    double rBp(void) const;
    double rBm(void) const;
    double xp(void) const;
    double xm(void) const;
    double yp(void) const;
    double ym(void) const;
    double gamma(void) const;
    double delb(void) const;

 private:
    void SetParamsGD(const std::vector<double>& par);
    void SetParamsXY(const std::vector<double>& par);
    bool CheckPars(const std::vector<double>& par, unsigned size) const;

    std::complex<double> m_zp;
    std::complex<double> m_zm;
    int m_flv;

    static const double radToDegrees;
    static const std::complex<double> imone;
};

#endif  // SRC_BTODK_DTOKSPIPI_MODEL_H_

