#include "btodk_dtokspipi_model.h"

#include <cmath>
#include <iostream>

using namespace std;

BtoDK_DtoKspipi_Model::BtoDK_DtoKspipi_Model(const double& gamma,const double& delb, const double& rb):
  KspipiModel(),
  m_flv(1)
{
  SetGD(gamma,delb,rb);
}

compld BtoDK_DtoKspipi_Model::Amp(const double& mAB, const double& mAC) const {
  const compld ap = KspipiModel::Amp(mAB,mAC);
  const compld an = KspipiModel::Amp(mAC,mAB);
  const compld z = m_flv>0 ? m_zp : m_zm;
//  cout << "Amp " << ap << " " << z << " " << an << endl;
  return ap + z*an;
}

void BtoDK_DtoKspipi_Model::SetFlv(const int flv) {m_flv = flv;}

void BtoDK_DtoKspipi_Model::SetGD(const double& gamma, const double& delta, const double& rB){
  m_zp = rB*exp(imone*(delta+gamma)/radToDegrees);
  m_zm = rB*exp(imone*(delta-gamma)/radToDegrees);
}

void BtoDK_DtoKspipi_Model::SetXY(const double& xp, const double& yp,const double& xm, const double& ym){
  m_zp = xp + imone*yp;
  m_zm = xm + imone*ym;
//  cout << "SetXY: " << m_zp << " " << m_zm << endl;
}

int    BtoDK_DtoKspipi_Model::Flv(void) const { return m_flv;}
double BtoDK_DtoKspipi_Model::rBp(void) const { return abs(m_zp);}
double BtoDK_DtoKspipi_Model::rBm(void) const { return abs(m_zm);}
double BtoDK_DtoKspipi_Model::xp( void) const { return real(m_zp);}
double BtoDK_DtoKspipi_Model::xm( void) const { return real(m_zm);}
double BtoDK_DtoKspipi_Model::yp( void) const { return imag(m_zp);}
double BtoDK_DtoKspipi_Model::ym( void) const { return imag(m_zm);}
double BtoDK_DtoKspipi_Model::gamma(void) const {
  return 0.5*arg(m_zp/m_zm)*radToDegrees;
}
double BtoDK_DtoKspipi_Model::delb(void)  const {
  return 0.5*arg(m_zp*m_zm)*radToDegrees;
}

void BtoDK_DtoKspipi_Model::SetParams(const vector<double>& par){
  switch(par.size()){
  case 3: SetParamsGD(par); break;
  case 4: SetParamsXY(par); break;
  default:
    cout << "BtoDK_DtoKspipi_Model::SetParams: wrong number of parameters " << par.size() << endl;
  }
}

int BtoDK_DtoKspipi_Model::SetParamsGD(const vector<double>& par){
  SetGD(par[0],par[1],par[2]);
  return 0;
}

int BtoDK_DtoKspipi_Model::SetParamsXY(const vector<double>& par){
  SetXY(par[0],par[1],par[2],par[3]);
  return 0;
}

