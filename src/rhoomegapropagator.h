/** Copyright 2017 Vitaly Vorobyev
 ** @file rhoomegapropagator.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_RHOOMEGAPROPAGATOR_H_
#define SRC_RHOOMEGAPROPAGATOR_H_

#include <complex>

#include "./gounarissakurai.h"
#include "./relbreitwigner.h"

/// \brief Class implementing rho(770)-omega(782) interfering propagator.
/// Taken from Eq.(13) in Phys. Rev. D92, 032002 (2015). Introduced for the
/// first time in R. R. Akhmetshin et al. (CMD-2 Collaboration),
/// Measurement of e+e- -> pi+ pi- cross section with CMD-2 around rho-meson,
/// Phys. Lett. B 527, 161 (2002).
class RhoOmegaPropagator : public GounarisSakurai {
 public:
    RhoOmegaPropagator(const double& a, const double& theta);

    std::complex<double> operator()(const double& s, const double& p) const;

 private:
    std::complex<double> m_omega_amp;
    RelBreitWigner* omega_prop;

    static const double m_PI_Mass;
    static const double m_rho770_Mass;
    static const double m_rho770_Width;
    static const double m_omega_Mass;
    static const double m_omega_Width;
    static const double m_p0_rho770;
    static const double m_p0_omega;
};

#endif  // SRC_RHOOMEGAPROPAGATOR_H_
