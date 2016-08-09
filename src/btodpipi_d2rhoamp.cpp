#include "btodpipi_d2rhoamp.h"

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

typedef BtoDpipi_D2rhoAmp  DrhoAmp;
typedef complex<double> compld;
typedef vector<double>  vectd;
typedef vector<compld>  vectcompld;

const compld ione(0,1);

cdouble mB_init  = 5.279;
cdouble mD_init  = 1.865;
cdouble mpi_init = 0.139568;

BtoDpipi_D2rhoAmp::BtoDpipi_D2rhoAmp():
  AbsDalitzModel(mB_init,mD_init,mpi_init,mpi_init),
  m_ffB(1.6), m_ffR(1.0),
  m_Ddst_mass(2.4677), m_Ddst_width(0.056002),
  m_rho_mass( 0.7756), m_rho_width( 0.15),
  m_Ddst_amp(1.,0.),
  m_rho_amp(0.428140*cos(1.43780),0.428140*sin(1.43780))
{
  SetResNames({"D2","rho"});
  SetAmpNames({"D2","rho"});
  SetABaxis("m(D#pi^{+}) (GeV^{2}/c^{4})");
  SetACaxis("m(D#pi^{-}) (GeV^{2}/c^{4})");
  SetBCaxis("m(#pi^{+}#pi^{-}) (GeV^{2}/c^{4})");
  init();
}

void DrhoAmp::init(void){
  SetCoefficients({m_Ddst_amp,m_rho_amp});
  SetAmpSignature({1,1});

  const double nsig = 0.5;
  vectd ledge = {
    pow(m_Ddst_mass-nsig*m_Ddst_width,2),
    pow(m_rho_mass - 0.5*m_rho_width, 2)
  };
  vectd redge = {
    pow(m_Ddst_mass+nsig*m_Ddst_width,2),
    pow(m_rho_mass + 0.5*m_rho_width, 2)
  };
  vecti types {1,3};
  SetResAreas(ledge,redge,types);
}

compld DrhoAmp::GetResVal(cdouble& mABsq,cdouble& mACsq, const int resnum) const {
  double pq  = 0;
  double hel = 0;
  cdouble mBCsq = GetmBCsq(mABsq,mACsq);
  if(resnum == 0) hel = CosHelAB(mABsq,mACsq,pq);// D pi resonance
  else            hel = CosHelBC(mBCsq,mABsq,pq);// pi pi resonance
  switch(resnum){
  case 0: return DdstAmp(mABsq,hel,pq); // D2
  case 1: return rhoAmp( mBCsq,hel,pq); // rho
  default: cout << "DrhoAmp: wrong resnum " << resnum << endl;
  }
  return 0;
}

void DrhoAmp::GetResVals(vectcd& resv, cdouble& mABsq,cdouble& mACsq) const {
  resv.clear();
  cdouble mBCsq = GetmBCsq(mABsq,mACsq);
  double pqAB;
  double pqBC;
  cdouble hAB = CosHelAB(mABsq,mACsq,pqAB);
  cdouble hBC = CosHelBC(mBCsq,mABsq,pqBC);
  resv.push_back(DdstAmp(mABsq,hAB,pqAB));// D** amplitude
  resv.push_back(rhoAmp( mBCsq,hBC,pqBC));// rho(770) amplitude
}

double DrhoAmp::Norm(void) const{ return NormWithCache();}

void DrhoAmp::SetParams(const vectd& pars){
  // parameters -> amplitudes and phases
  if(pars.size() != 2*AmpNum()){
    cout << "Wrong pars size: " << pars.size() << " (" << 2*AmpNum() << " expected)" << endl;
    return;
  }
  vectcd coeffs;
  for(unsigned i=0; i<AmpNum(); i++){
    cdouble amp = pars[2*i];
    cdouble pha = pars[2*i+1];
    coeffs.push_back(amp*compld(cos(pha),sin(pha)));
  }
  SetCoefficients(coeffs);
  return;
}

double DrhoAmp::ang1(cdouble& hel,cdouble& pq) const { return pq*hel;}
double DrhoAmp::ang2(cdouble& hel,cdouble& pq) const { return pow(pq,2)*(pow(hel,2)-1./3.);}

compld DrhoAmp::DdstAmp(cdouble& mABsq,cdouble& helAB,cdouble& pqAB) const{
  return ang2(helAB,pqAB)*AmpABspin2(mABsq,m_Ddst_mass,m_Ddst_width);
}

compld DrhoAmp::rhoAmp(cdouble& mBCsq,cdouble& helBC,cdouble& pqBC) const{
    return ang1(helBC,pqBC)*AmpBCspin1(mBCsq,m_rho_mass,m_rho_width);
}

double DrhoAmp::BWUnit1(cdouble& p, cdouble& scale) const{
  return 1. + pow(p*scale,2);
}
double DrhoAmp::BWUnit2(cdouble& p, cdouble& scale) const{
  cdouble psc = p*scale;
  return 9. + 3.*pow(psc,2) + pow(psc,4);
}

double DrhoAmp::BlattWeisskopf(cdouble& p0,cdouble& p, cint spin, cdouble& scale) const{
  switch(spin){
  case 0: return 1.;
  case 1: return sqrt(BWUnit1(p0,scale)/BWUnit1(p,scale));
  case 2: return sqrt(BWUnit2(p0,scale)/BWUnit2(p,scale));
  }
  return 1.;
}

compld DrhoAmp::AmpABspin2(cdouble& mABsq,cdouble& mR,cdouble& gR) const{
  cdouble mRsq = pow(mR,2);
  // pion momentum in resonance frame
  cdouble ppi0 = pB_AB(mRsq);
  cdouble ppi  = pB_AB(mABsq);
  // resonance momentum in B0 frame
  cdouble pR0  = pAB_M(mRsq);
  cdouble pR   = pAB_M(mABsq);
  // B0 formfactor
  cdouble ffB  = BlattWeisskopf(pR0,pR,m_ffB,2);
  // resonance formfactor
  cdouble ffR  = BlattWeisskopf(ppi0,ppi,m_ffR,2);
  // Breit-Wigner
  cdouble ar = mABsq-mRsq;
  cdouble gr = (gR*mRsq/sqrt(mABsq))*pow(ppi/ppi0,5)*pow(ffB,2);
  return ffB*ffR/compld(ar,gr);
}

compld DrhoAmp::AmpBCspin1(cdouble& mBCsq,cdouble& mR,cdouble& gR) const{
  cdouble mRsq = pow(mR,2);
  cdouble ppi0 = pB_BC(mRsq);
  cdouble ppi  = pB_BC(mBCsq);
  if(ppi == 0){
     cout << "AmpBCspin1: E(pi) < m(pi) " << eB_BC(mBCsq) << " " << mB() << " " << mBCsq << endl;
     return 0;
  }
  double hwd0,hwd;
  cdouble ar  = mBCsq-mRsq;
  cdouble gr  = (gR*mRsq/sqrt(mBCsq))*pow(ppi/ppi0,3);
  cdouble hw  = hwrho(mBCsq,hwd);
  cdouble hw0 = hwrho(mRsq,hwd0);
  cdouble dm  = ((hw-hw0)*pow(ppi,2)-ar*hwd0*pow(ppi0,2))*mRsq*gR/pow(ppi0,3);
  return 1./compld(ar-dm,gr);
}

double DrhoAmp::hwrho(cdouble& s, double& hwd) const{
  cdouble y = sqrt(1.-0.07795/s);
  cdouble w = log((1.+y)/(1-y));
  hwd = (0.5*w*(1.-pow(y,2))/y+1.)/s*M_PI_2;
  return w*y*M_PI_2;
}
