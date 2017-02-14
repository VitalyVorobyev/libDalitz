/** Copyright 2017 Vitaly Vorobyev
 ** @file kspipimodel.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/kspipimodel.h"

#include <cmath>
#include <complex>

#include "./blattweisskopf.h"
#include "./absvarwidth.h"

const int AB = DalitzResonance::AB;
const int BC = DalitzResonance::BC;
const int AC = DalitzResonance::AC;

const double KspipiModel::dtr = M_PI / 180.;
const double KspipiModel::m_D0_Mass = 1.865;
const double KspipiModel::m_Ks0_Mass = 0.497611;
const double KspipiModel::m_PI_Mass = 0.13957018;

KspipiModel::KspipiModel(void) :
    KspipiModel(m_D0_Mass, m_Ks0_Mass, m_PI_Mass) {}

KspipiModel::KspipiModel(const double &md, const double &mks,
                         const double& mpi) :
    SymDalitzModel(md, mks, mpi, -M_PI/8., 15.*M_PI/8.) {
    SetABaxis("m_{+}^{2}, GeV^{2}/c^{4}");
    SetACaxis("m_{-}^{2}, GeV^{2}/c^{4}");
    SetBCaxis("m_{#pi#pi}^{2}, GeV^{2}/c^{4}");
    // ** A. Poluektov et al. Phys. Rev. D 81, 112002 – 16 June 2010 **
    BlattWeisskopf::m_r_meson     = 5.0;
    BlattWeisskopf::m_r_resonance = 1.5;
    ResDecayAngularDistribution::m_use_mRsq = true;

    // CF //
    AddRes(new DalitzResonance("K*(892)", ResPropType::RBW, VarWType::BW,
               this, AB, 0.8937, 0.0484, 1, 1.638, 133.2*dtr));
    AddRes(new DalitzResonance("K0*(1430)", ResPropType::RBW,
               this, AB, 1.4120, 0.2940, 0, 2.210, 358.9*dtr));
    AddRes(new DalitzResonance("K2*(1430)", ResPropType::RBW,
               this, AB, 1.4256, 0.0985, 2, 0.890, 314.8*dtr));
    AddRes(new DalitzResonance("K*(1680)", ResPropType::RBW,
               this, AB, 1.7170, 0.3220, 1, 0.880, 82.0*dtr));
    AddRes(new DalitzResonance("K*(1410)", ResPropType::RBW,
               this, AB, 1.4140, 0.2320, 1, 0.650, 120.0*dtr));
    // DCS //
    AddRes(new DalitzResonance("K*(892) DCS", ResPropType::RBW,
               this, AC, .89166, 0.0508, 1, 0.149, 325.4*dtr));
    AddRes(new DalitzResonance("K0*(1430) DCS", ResPropType::RBW,
               this, AC, 1.4120, 0.2940, 0, 0.360, 87.0*dtr));
    AddRes(new DalitzResonance("K2*(1430) DCS", ResPropType::RBW,
               this, AC, 1.4256, 0.0985, 2, 0.230, 275.0*dtr));
    AddRes(new DalitzResonance("K*(1680) DCS", ResPropType::RBW,
               this, AC, 1.7170, 0.3220, 1, 2.100, 130.0*dtr));
    AddRes(new DalitzResonance("K*(1410) DCS", ResPropType::RBW,
               this, AC, 1.4140, 0.2320, 1, 0.420, 253.0*dtr));
    // CP
    AddRes(new DalitzResonance("rho(770)", ResPropType::RBW,
               this, BC, 0.7717, 0.1490, 1, 1.000, 0.000));
    AddRes(new DalitzResonance("omega(782)", ResPropType::RBW,
               this, BC, .78265, .00849, 1, .0343, 112.0*dtr));
    AddRes(new DalitzResonance("f0(980)", ResPropType::RBW,
               this, BC, 0.9770, 0.0500, 0, 0.385, 207.3*dtr));
    AddRes(new DalitzResonance("f0(1370)", ResPropType::RBW,
               this, BC, 1.3100, 0.2720, 0, 1.250, 69.0*dtr));
    AddRes(new DalitzResonance("f2(1270)", ResPropType::RBW,
               this, BC, 1.2754, 0.1851, 2, 1.440, 342.9*dtr));
    AddRes(new DalitzResonance("rho(1450)", ResPropType::RBW,
               this, BC, 1.4650, 0.4000, 1, 0.490, 64.0*dtr));
    AddRes(new DalitzResonance("sigma1", ResPropType::RBW,
               this, BC, 0.5220, 0.4530, 0, 1.560, 214.0*dtr));
    AddRes(new DalitzResonance("sigma2", ResPropType::RBW,
               this, BC, 1.0330, 0.0880, 0, 0.200, 212.0*dtr));
    // NR
    AddRes(new DalitzResonance("NR", ResPropType::NR,
               this, BC, 0, std::complex<double>(-2.537, 0.923)));
}
