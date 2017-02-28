/** Copyright 2017 Vitaly Vorobyev
 ** @file btodk_dtokspipi_model.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/btodk_dtokspipi_model.h"

#include <cmath>
#include <iostream>

typedef BtoDK_DtoKspipi_Model BDK;
typedef std::complex<double> compld;

using std::cout;
using std::endl;

using std::abs;
using std::real;
using std::imag;
using std::exp;

const double BDK::radToDegrees =  180./M_PI;
const compld BDK::imone = compld(0, 1);

BDK::BtoDK_DtoKspipi_Model(
        const double& gamma, const double& delb, const double& rb) :
    AbsDalitzModel(0,0,0,0),
    KspipiModel(),
    m_flv(1) {
    SetGD(gamma, delb, rb);
}

compld BDK::Amp(const double& mAB, const double& mAC) const {
    const compld ap = KspipiModel::Amp(mAB, mAC);
    const compld an = KspipiModel::Amp(mAC, mAB);
    const compld z = m_flv > 0 ? m_zp : m_zm;
    return ap + z*an;
}

void BDK::SetFlv(const int flv) {m_flv = flv;}

void BDK::SetGD(const double& gamma, const double& delta,
                                  const double& rB) {
    m_zp = rB*exp(imone*(delta+gamma)/radToDegrees);
    m_zm = rB*exp(imone*(delta-gamma)/radToDegrees);
}

void BDK::SetXY(const double& xp, const double& yp,
                                  const double& xm, const double& ym) {
    m_zp = xp + imone*yp;
    m_zm = xm + imone*ym;
}

int BDK::Flv(void) const { return m_flv;}
double BDK::rBp(void) const { return abs(m_zp);}
double BDK::rBm(void) const { return abs(m_zm);}
double BDK::xp( void) const { return real(m_zp);}
double BDK::xm( void) const { return real(m_zm);}
double BDK::yp( void) const { return imag(m_zp);}
double BDK::ym( void) const { return imag(m_zm);}

double BDK::gamma(void) const {
  return 0.5*arg(m_zp/m_zm)*radToDegrees;
}

double BDK::delb(void)  const {
  return 0.5*arg(m_zp*m_zm)*radToDegrees;
}

void BDK::SetParams(const std::vector<double>& par) {
    switch (par.size()) {
    case 3: SetParamsGD(par); break;
    case 4: SetParamsXY(par); break;
    default:
        cout << "BtoDK_DtoKspipi_Model::SetParams: wrong number of "
             << "parameters " << par.size() << endl;
    }
}

void BDK::SetParamsGD(const std::vector<double>& par) {
    if (!CheckPars(par, 3)) return;
    SetGD(par[0], par[1], par[2]);
}

void BDK::SetParamsXY(const std::vector<double>& par) {
    if (!CheckPars(par, 4)) return;
    SetXY(par[0], par[1], par[2], par[3]);
}

bool BDK::CheckPars(const std::vector<double>& par, unsigned size) const {
    if (par.size() != size) {
        cout << "BDK::SetParamsGD: wrong par vector size: "
             << par.size() << endl;
        cout << size << " expected" << endl;
        return false;
    }
    return true;
}

