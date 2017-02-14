/** Copyright 2017 Vitaly Vorobyev
 ** @file btodpipi_d2rhoamp.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/btodpipi_d2rhoamp.h"

#include <iostream>
#include <cmath>
#include <fstream>

typedef BtoDpipi_D2rhoAmp DrhoAmp;
typedef std::complex<double> compld;
typedef std::vector<double> vectd;
typedef std::vector<int> vecti;
typedef std::vector<compld> vectcd;

using std::cout;
using std::endl;
using std::sqrt;
using std::pow;
using std::log;
using std::sin;
using std::cos;
using std::fabs;
//using std::M_PI_2;

const double DrhoAmp::mB0 = 5.279;
const double DrhoAmp::mD0 = 1.865;
const double DrhoAmp::mpi = 0.139568;
const double DrhoAmp::rB = 1.6;
const double DrhoAmp::rR = 1.0;
const double DrhoAmp::Ddst_mass = 2.4677;
const double DrhoAmp::Ddst_width = 0.056002;
const double DrhoAmp::rho_mass = 0.7756;
const double DrhoAmp::rho_width = 0.15;
const compld DrhoAmp::Ddst_amp = compld(1., 0.);
const compld DrhoAmp::rho_amp = compld(0.428140*cos(1.43780), 0.428140*sin(1.43780));

BtoDpipi_D2rhoAmp::BtoDpipi_D2rhoAmp() :
    AbsDalitzModel(mB0, mD0, mpi, mpi) {
    SetResNames({"D2", "rho"});
    SetAmpNames({"D2", "rho"});
    SetABaxis("m(D#pi^{+}) (GeV^{2}/c^{4})");
    SetACaxis("m(D#pi^{-}) (GeV^{2}/c^{4})");
    SetBCaxis("m(#pi^{+}#pi^{-}) (GeV^{2}/c^{4})");
    init();
}

void DrhoAmp::init(void) {
    SetCoefficients({Ddst_amp, rho_amp});
    SetAmpSignature({1, 1});
    SetParState({false, true, true, true});

    const double nsig = 0.5;
    vectd ledge = {
        pow(Ddst_mass-nsig*Ddst_width, 2),
        pow(rho_mass - 0.5*rho_width, 2)
    };
    vectd redge = {
        pow(Ddst_mass+nsig*Ddst_width, 2),
        pow(rho_mass + 0.5*rho_width, 2)
    };
    vecti types {1, 3};
    SetResAreas(ledge, redge, types);
}

compld DrhoAmp::GetResVal(const double& mABsq, const double& mACsq,
                          const int resnum) const {
    double pq  = 0;
    double hel = 0;
    const double mBCsq = GetmBCsq(mABsq, mACsq);
    if (resnum == 0) hel = CosHelAB(mABsq, mACsq, &pq);  // D pi resonance
    else             hel = CosHelBC(mBCsq, mABsq, &pq);  // pi pi resonance
    switch (resnum) {
    case 0: return DdstAmp(mABsq, hel, pq);  // D2
    case 1: return rhoAmp(mBCsq, hel, pq);  // rho
    default: cout << "DrhoAmp: wrong resnum " << resnum << endl;
    }
    return 0;
}

void DrhoAmp::GetResVals(vectcd* resv, const double& mABsq,
                         const double& mACsq) const {
    resv->clear();
    const double mBCsq = GetmBCsq(mABsq, mACsq);
    double pqAB;
    double pqBC;
    const double hAB = CosHelAB(mABsq, mACsq, &pqAB);
    const double hBC = CosHelBC(mBCsq, mABsq, &pqBC);
    resv->push_back(DdstAmp(mABsq, hAB, pqAB));  // D** amplitude
    resv->push_back(rhoAmp(mBCsq, hBC, pqBC));  // rho(770) amplitude
}

double DrhoAmp::Norm(void) const { return NormWithCache();}

void DrhoAmp::SetParams(const vectd& pars) {
    // parameters -> amplitudes and phases
    if (pars.size() != 2*AmpNum()) {
        cout << "Wrong pars size: " << pars.size()
             << " (" << 2*AmpNum() << " expected)" << endl;
        return;
    }
    vectcd coeffs;
    for (unsigned i=0; i < AmpNum(); i++) {
        const double amp = pars[2*i];
        const double pha = pars[2*i+1];
        coeffs.push_back(amp*compld(cos(pha), sin(pha)));
    }
    SetCoefficients(coeffs);
}

double DrhoAmp::ang1(const double& hel, const double& pq) const {
    return pq*hel;
}
double DrhoAmp::ang2(const double& hel, const double& pq) const {
    return pow(pq, 2)*(pow(hel, 2)-1./3.);
}

compld DrhoAmp::DdstAmp(const double& mABsq, const double& helAB,
                        const double& pqAB) const {
  return ang2(helAB, pqAB)*AmpABspin2(mABsq, Ddst_mass, Ddst_width);
}

compld DrhoAmp::rhoAmp(const double& mBCsq, const double& helBC,
                       const double& pqBC) const {
    return ang1(helBC, pqBC)*AmpBCspin1(mBCsq, rho_mass, rho_width);
}

double DrhoAmp::BWUnit1(const double& p, const double& scale) const {
  return 1. + pow(p*scale, 2);
}

double DrhoAmp::BWUnit2(const double& p, const double& scale) const {
  const double psc = p*scale;
  return 9. + 3.*pow(psc, 2) + pow(psc, 4);
}

double DrhoAmp::BlattWeisskopf(const double& p0, const double& p,
                               const int spin, const double& scale) const {
    switch (spin) {
    case 0: return 1.;
    case 1: return sqrt(BWUnit1(p0, scale)/BWUnit1(p, scale));
    case 2: return sqrt(BWUnit2(p0, scale)/BWUnit2(p, scale));
    }
    return 1.;
}

compld DrhoAmp::AmpABspin2(const double& mABsq, const double& mR,
                           const double& gR) const {
    const double mRsq = pow(mR, 2);
    // pion momentum in resonance frame
    const double ppi0 = pB_AB(mRsq);
    const double ppi  = pB_AB(mABsq);
    // resonance momentum in B0 frame
    const double pR0  = pAB_M(mRsq);
    const double pR   = pAB_M(mABsq);
    // B0 formfactor
    const double ffB  = BlattWeisskopf(pR0, pR, rB, 2);
    // resonance formfactor
    const double ffR  = BlattWeisskopf(ppi0, ppi, rR, 2);
    // Breit-Wigner
    const double ar = mABsq-mRsq;
    const double gr = (gR*mRsq/sqrt(mABsq))*pow(ppi/ppi0, 5)*pow(ffB, 2);
    return ffB*ffR/compld(ar, gr);
}

compld DrhoAmp::AmpBCspin1(const double& mBCsq, const double& mR,
                           const double& gR) const {
    const double mRsq = pow(mR, 2);
    const double ppi0 = pB_BC(mRsq);
    const double ppi  = pB_BC(mBCsq);
    if (fabs(ppi) < 0.0001) {
        cout << "AmpBCspin1: E(pi) < m(pi) " << eB_BC(mBCsq)
             << " " << mB() << " " << mBCsq << endl;
        return 0;
    }
    double hwd0, hwd;
    const double ar  = mBCsq-mRsq;
    const double gr  = (gR*mRsq/sqrt(mBCsq))*pow(ppi/ppi0, 3);
    const double hw  = hwrho(mBCsq, &hwd);
    const double hw0 = hwrho(mRsq, &hwd0);
    const double dm  = ((hw-hw0)*pow(ppi, 2)-ar*hwd0*pow(ppi0, 2))*
            mRsq*gR/pow(ppi0, 3);
    return 1./compld(ar-dm, gr);
}

double DrhoAmp::hwrho(const double& s, double *hwd) const {
    const double y = sqrt(1.-0.07795/s);
    const double w = log((1.+y)/(1-y));
    *hwd = (0.5*w*(1.-pow(y, 2))/y+1.)/s*M_PI_2;
    return w*y*M_PI_2;
}
