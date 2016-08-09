#ifndef CONSTS_H
#define CONSTS_H

#include <complex>
#include <cmath>

//const double EvtConst::pi           =  3.141592653589793238;
//const double EvtConst::twoPi        =  2*pi;
const double radToDegrees =  180./M_PI;
//const double EvtConst::c            =  2.99792458E11;       // mm/sec

const double m_PI_Mass  = 0.13957018;
const double m_PI0_Mass = 0.1349767;
const double m_Ks0_Mass = 0.497614;
const double m_D0_Mass  = 1.86484;
const double m_B0_Mass  = 5.27958;
const double m_K_Mass   = 0.493667;
const double m_eta_Mass = 0.547862;

// Phys. Rev. D 86, 032013 (2012)
const double m_rho770_Mass   = 0.77502;// +- 0.00035
const double m_rho770_Width  = 0.14959;// +- 0.00067

const double m_omega_Mass    = 0.78191;// +- 0.00024
const double m_omega_Width   = 0.00813;// +- 0.00045

const double m_rho1450_Mass  = 1.493;// +- 0.015
const double m_rho1450_Width = 0.427;// +- 0.031

const double m_rho1700_Mass  = 1.861;// +- 0.017
const double m_rho1700_Width = 0.316;// +- 0.026
//

typedef std::complex<double> compld;
const compld imone = compld(0,1);

#endif // CONSTS_H
