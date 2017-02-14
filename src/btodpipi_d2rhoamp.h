/** Copyright 2017 Vitaly Vorobyev
 ** @file btodpipi_d2rhoamp.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_BTODPIPI_D2RHOAMP_H_
#define SRC_BTODPIPI_D2RHOAMP_H_

#include <vector>
#include <complex>

#include "./absdalitzmodel.h"

///
/// \brief The BtoDpipi_D2rhoAmp class
///
class BtoDpipi_D2rhoAmp : public AbsDalitzModel {
 public:
    BtoDpipi_D2rhoAmp();
    /// Vector of Res_amp (rho, omega and rho' combined)
    void GetResVals(std::vector<std::complex<double>>* resv,
                    const double& mABsq, const double& mACsq) const;
    /// Complex amplitude of a resonance
    std::complex<double> GetResVal(const double& mABsq,
                                   const double& mACsq, const int resnum) const;
    /// Set parameters
    void SetParams(const std::vector<double>& pars);
    /// Get normalization
    double Norm(void) const;

 private:
    void init(void);

    std::complex<double> AmpABspin2(const double& mABsq,
                                    const double& mR, const double& gR) const;
    std::complex<double> AmpBCspin1(const double& mBCsq,
                                    const double& mR, const double& gR) const;

    std::complex<double> DdstAmp(const double& mABsq,
                                 const double& helAB, const double& pqAB) const;
    std::complex<double> rhoAmp(const double& mBCsq,
                                const double& helBC, const double& pqBC) const;

    double hwrho(const double& s, double* hwd) const;

    double ang1(const double& hel, const double& pq) const;
    double ang2(const double& hel, const double& pq) const;

    double BWUnit1(const double& p, const double& scale) const;
    double BWUnit2(const double& p, const double& scale) const;
    double BlattWeisskopf(const double& p0, const double& p,
                          const int spin, const double& scale) const;

    static const double mB0;
    static const double mD0;
    static const double mpi;

    static const double rB;
    static const double rR;

    static const double Ddst_mass;
    static const double Ddst_width;
    static const double rho_mass;
    static const double rho_width;

    static const std::complex<double> Ddst_amp;
    static const std::complex<double> rho_amp;
};

#endif  // SRC_BTODPIPI_D2RHOAMP_H_
