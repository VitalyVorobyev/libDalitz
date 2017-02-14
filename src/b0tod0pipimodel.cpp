/** Copyright 2017 Vitaly Vorobyev
 ** @file b0tod0pipimodel.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/b0tod0pipimodel.h"

#include <cmath>

#include "./blattweisskopf.h"
#include "./dalitzvetodst2010.h"

const double B0toD0pipiModel::m_B0_Mass = 5.279;
const double B0toD0pipiModel::m_D0_Mass = 1.865;
const double B0toD0pipiModel::m_PI_Mass = 0.139568;
const double B0toD0pipiModel::dtr = M_PI / 180.;

const int AB = DalitzResonance::AB;
const int BC = DalitzResonance::BC;

const int B0toD0pipiModel::Belle = 0;
const int B0toD0pipiModel::LHCb = 1;

B0toD0pipiModel::B0toD0pipiModel(const int type) :
    B0toD0pipiModel(m_B0_Mass, m_D0_Mass, m_PI_Mass, type) {}

B0toD0pipiModel::B0toD0pipiModel(const double& mB, const double& mD,
                                 const double& mpi, const int type) :
    SymDalitzModel(mB, mD, mpi, -M_PI/8., 15.*M_PI/8.), m_type(type) {
    SetABaxis("m^{2}(D^{0}#pi^{+}), GeV^{2}/c^{4}");
    SetACaxis("m^{2}(D^{0}#pi^{-}), GeV^{2}/c^{4}");
    SetBCaxis("m^{2}(#pi^{+}#pi^{-}), GeV^{2}/c^{4}");
    if (m_type == Belle) InitBelleModel();
    if (m_type == LHCb) InitLHCbModel();
}

void B0toD0pipiModel::InitLHCbModel(void) {
    BlattWeisskopf::m_r_meson = 1.6;
    BlattWeisskopf::m_r_resonance = 1.6;
    //                             mass     width    J   amp     phase
    AddRes(new DalitzResonance("Dv*(2010)", ResPropType::VDst, this,
                     AB,  0.95,    0.51,       18.80, 266.7*dtr));
    AddRes(new DalitzResonance("D*0(2400)", ResPropType::RBW, this,
                     AB,  2.349,   0.217,   0, 12.10,  83.6*dtr));
    AddRes(new DalitzResonance("D*2(2460)", ResPropType::RBW, this,
                     AB,  2.4686,  0.0473,  2,  1.310, 262.9*dtr));
    AddRes(new DalitzResonance("D*J(2760)", ResPropType::RBW, this,
                     AB,  2.798,   0.105,   3,  0.053,  91.1*dtr));
    AddRes(new DalitzResonance("rho(770)",  ResPropType::RhoOmega, this,
                     BC,  0.30,    176.8*dtr,   1.000,   0.0*dtr));
//    AddRes(new DalitzResonance("rho(770)",  ResPropType::GS, this,
//                     BC,  0.7700,  0.1300,  1,  1.000,   0.0*dtr));
//    AddRes(new DalitzResonance("rho(770)",  ResPropType::RBW, this,
//                     BC,  0.77502, 0.14959, 1,  1.000,   0.0*dtr));
    AddRes(new DalitzResonance("rho(1450)", ResPropType::GS, this,
                     BC,  1.4930,  0.4270,  1,  0.230, 149.0*dtr));
    AddRes(new DalitzResonance("rho(1700)", ResPropType::GS, this,
                     BC,  1.8610,  0.3160,  1,  0.078, 103.5*dtr));
    AddRes(new DalitzResonance("f2(1270)",  ResPropType::RBW, this,
                     BC,  1.2751,  0.1851,  2,  0.072, 158.1*dtr));
    AddRes(new DalitzResonance("f0(500)",   ResPropType::Bugg, this,
                     BC,                       18.70,   38.4*dtr));
    AddRes(new DalitzResonance("f0(980)",   ResPropType::Flatte, this,
                     BC,  0.9399,               2.620, 138.9*dtr));
    AddRes(new DalitzResonance("f0(2020)",  ResPropType::RBW, this,
                     BC,  1.992,   0.442,   0,  4.410, 258.5*dtr));
    AddRes(new DalitzResonance("NR",        ResPropType::NR, this,
                     BC, -0.363,                3.430,  77.1*dtr));
}

void B0toD0pipiModel::InitBelleModel(void) {
    BlattWeisskopf::m_r_meson     = 1.6;
    BlattWeisskopf::m_r_resonance = 1.6;
    const double omegaAmp = 2.0*0.027/(0.7756*0.7756);
    const double omegaPha = 1.99+2.25;
  //                     mass  width      J amp   phase
    AddRes(new DalitzResonance("D2*", ResPropType::RBW, this,
                     AB, 2.46570, 0.04960, 2, 1.00,  0.00));
    AddRes(new DalitzResonance("D0*", ResPropType::RBW, this,
                     AB, 2.30800, 0.27611, 0, 40.6, -3.00));
    AddRes(new DalitzResonance("Dv*", ResPropType::VDst2, this,
                     AB,                  271., -2.62));
    AddRes(new DalitzResonance("Dv*", ResPropType::RBW, this,
                     AB, 2.0103, 83.4e-6, 1, 0.0000, -2.62));
    AddRes(new DalitzResonance("rho(770)", ResPropType::GS, this,
                     BC, 0.77560, 0.14400, 1, 2.0, 2.25));
    AddRes(new DalitzResonance("rho(770)", ResPropType::RBW, this,
                     BC, 0.7756, 0.14400, 1, 2.0, 2.25));
    AddRes(new DalitzResonance("omega",    ResPropType::RBW, this,
                     BC, 0.78265, 0.00890, 1, omegaAmp, omegaPha));
    AddRes(new DalitzResonance("f2", ResPropType::RBW, this,
                     BC, 1.27500, 0.18500, 2, 0.10, 2.97));
    AddRes(new DalitzResonance("f0(600)",  ResPropType::RBW, this,
                     BC, 0.51300, 0.33500, 0, 15.7, -0.44));
    AddRes(new DalitzResonance("f0(980)",  ResPropType::RBW, this,
                     BC, 0.97800, 0.04400, 0, 3.08, -2.48));
    AddRes(new DalitzResonance("f0(1370)", ResPropType::RBW, this,
                     BC, 1.43400, 0.17300, 0, 10.2, -1.52));
    // Cut off 3 MeV around D*(2010) mass
    AddVeto(new DalitzVetoDst2010(0.003));
}

// START MIGRAD MINIMIZATION.  STRATEGY 1.  CONVERGENCE WHEN EDM .LT. 0.10E-03

// FCN=   112567.1     FROM MIGRAD    STATUS=INITIATE     254 CALLS  256 TOTAL
//                     EDM= unknown      STRATEGY= 1      NO ERROR MATRIX

//  EXT PARAMETER               CURRENT GUESS       STEP         FIRST
//  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE
//   1      a0       0.50000E-01   0.10000E-01    0.0000       -8044.1
//   2      a1       0.90000E-01   0.10000E-01    0.0000       -3438.1
//   3      f0       -3.6039       0.10000E-01    0.0000        336.39
//   4      f1        3.1074       0.10000E-01    0.0000        150.22
//   5      am2       2.4689       0.10000E-02    0.0000        8.1039
//   6      gm2      0.44170E-01   0.10000E-02    0.0000       -121.89
//   7      am0       2.3080       constant
//   8      gm0      0.27611       constant
//   9      evs       2000.0        1.0000        0.0000       -52.130
//  10      ab        0.0000       constant
//  11      fb        1.0101       0.10000        0.0000        3.1368
//  12      mr       0.77560       constant
//  13      gr       0.15000       constant
//  14      ar       0.29000       0.10000E-01    0.0000       -985.64
//  15      fr        1.7400       0.10000        0.0000        165.22
//  16      mf0      0.51300       constant
//  17      gf0      0.33500       constant
//  18      af0      0.64763E-01   0.10000E-01    0.0000       -252.33
//  19      ff0      -1.6046       0.10000        0.0000        255.28
//  20      mf2       1.2750       constant
//  21      gf2      0.18500       constant
//  22      af2      0.10000E+00   0.10000E-01    0.0000       -227.69
//  23      ff2       2.3400       0.10000        0.0000       -298.75
//  24      a3        1.9657       0.10000E-01   0.40463E-03   -913.69
//  25      f3       0.75967       0.10000        0.0000       -3.2656
//  26      fdf       1.6000       constant
//  27      fbf       1.0000       constant
//  28      a7        0.0000       constant
//  29      f7        0.0000       constant
//  30      amom     0.27200E-01   constant

// COVARIANCE MATRIX CALCULATED SUCCESSFULLY

// FCN=   112390.9     FROM HESSE     STATUS=OK       213 CALLS     1368 TOTAL
//                     EDM=  0.37E-03    STRATEGY= 1      ERROR MATRIX ACCURATE

//  EXT PARAMETER                                   STEP         FIRST
//  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE
//   1      a0       0.35749E-01   0.20165E-01   0.19461E-03   -5.9093
//   2      a1       0.77469E-01   0.10459E-01   0.14769E-03   -3.0071
//   3      f0       -3.8693       0.14167       0.30859E-02   0.34416
//   4      f1        2.7533       0.17147       0.32832E-02  -0.57016
//   5      am2       2.4677       0.19306E-02   0.42549E-01   0.29629E-01
//   6      gm2      0.56000E-01   0.52668E-02   0.20236E-01   0.75391E-01
//   7      am0       2.3080       constant
//   8      gm0      0.27611       constant
//   9      evs       2149.5        69.452       0.85965E-02   0.22247E-01
//  10      ab        0.0000       constant
//  11      fb        1.0150       0.38529       0.16890E-01   0.10383
//  12      mr       0.77560       constant
//  13      gr       0.15000       constant
//  14      ar       0.42748       0.23901E-01   0.33832E-02  -0.30510
//  15      fr        1.4390       0.11228       0.20965E-02   0.46707E-01
//  16      mf0      0.51300       constant
//  17      gf0      0.33500       constant
//  18      af0      0.55581E-01   0.10506E-01   0.56636E-02   0.11173
//  19      ff0      -1.7218       0.14463       0.28293E-02   0.24296
//  20      mf2       1.2750       constant
//  21      gf2      0.18500       constant
//  22      af2      0.81378E-01   0.99788E-02   0.44101E-02   0.27231
//  23      ff2       2.7154       0.14684       0.32218E-02  -0.16524
//  24      a3        1.9663       0.58796E-03   0.73267E-03    2.8269
//  25      f3       0.80074       0.34972       0.11907E-01  -0.15289
//  26      fdf       1.6000       constant
//  27      fbf       1.0000       constant
//  28      a7        0.0000       constant
//  29      f7        0.0000       constant
//  30      amom     0.27200E-01   constant

// br 1:   0. #   12.8338789
// br 2:   0.0343047469 #   5.7399365
// br 3:   0.100986009 #   11.1204642
// br 4:   0.410044718 #   1.51463363
// br 5:   0.0531335029 #   0.279036065
// br 6:   0.0694788905 #   0.339540311
// br 7:   0. #   0.
// br 8:   0.00207668497 #   0.00728203706
