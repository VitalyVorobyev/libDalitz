#include "kuzbtodpipiamp.h"

#include <iostream>
#include <cmath>
#include <fstream>

//#include<boost/range/numeric.hpp>

using namespace std;

typedef KuzBtoDpipiAmp  KAmp;
typedef complex<double> compld;
typedef vector<double>  vectd;
typedef vector<compld>  vectcompld;
typedef cdouble    cdouble;

const compld ione(0,1);

cdouble KAmp::m_ffr   = 1.6;
cdouble KAmp::m_ffdel = 1.0;

KAmp::KuzBtoDpipiAmp():
  AbsDalitzModel(m_fm,m_dstm,m_pm,m_pm)
{
  SetResNames({"D*","D2","D0","rho","omega","rho'","f2","f0(550)","f(980)","f(1370)"});
  SetAmpNames({"D*","D2","D0","rho","f2","f0(550)","f(980)","f(1370)"});
  SetABaxis("m(D#pi^{+}) (GeV^{2}/c^{4})");
  SetACaxis("m(D#pi^{-}) (GeV^{2}/c^{4})");
  SetBCaxis("m(#pi^{+}#pi^{-}) (GeV^{2}/c^{4})");
  init();
  PrintSummary();
}

compld KAmp::DstarAmp(  cdouble& qp,cdouble& h2,cdouble& pq2) { return ang1(h2,pq2)*amrs1(  qp);}
compld KAmp::DdstarAmp( cdouble& qp,cdouble& h2,cdouble& pq2) { return ang2(h2,pq2)*amrs2(  qp,m_ams2,m_gms2);}
compld KAmp::DzeroAmp(  cdouble& qp)                          { return              amrs0(  qp,m_ams0,m_gms0);}
compld KAmp::rhoAmp(    cdouble& q3,cdouble& h3,cdouble& pq3) { return ang1(h3,pq3)*amrs1h( q3,m_amrh,m_gmrh);}
compld KAmp::omegaAmp(  cdouble& q3,cdouble& h3,cdouble& pq3) { return ang1(h3,pq3)*amrs1hm(q3)*m_amom*exp(1.99*ione)*(-1.);}
compld KAmp::rhopAmp(   cdouble& q3,cdouble& h3,cdouble& pq3) { return ang1(h3,pq3)*amrs1h( q3,m_amrh1,m_gmrh1)*(-0.248);}
compld KAmp::f2Amp(     cdouble& q3,cdouble& h3,cdouble& pq3) { return ang2(h3,pq3)*amrs2h( q3,m_amf2,m_gmf2);}
compld KAmp::f0Amp(     cdouble& q3)                          { return              amrs0h( q3,m_amf0,m_gmf0);}
compld KAmp::f980Amp(   cdouble& q3)                          { return              amrs0h( q3,m_amf098,m_gmf098);}
compld KAmp::f1370Amp(  cdouble& q3)                          { return              amrs0h( q3,m_am3,m_gm3);}
compld KAmp::FullRhoAmp(cdouble& q3,cdouble& h3,cdouble& pq3){
  return ang1(h3,pq3)*(amrs1h( q3,m_amrh,m_gmrh)
                      +amrs1hm(q3)*m_amom*exp(1.99*ione)*(-1.)
                      +amrs1h( q3,m_amrh1,m_gmrh1)*(-0.248));
}

compld KAmp::GetResVal(cdouble& mABsq,cdouble& mACsq, const int resnum) const{
  double pq  = 0;
  double hel = 0;
  cdouble mBCsq = GetmBCsq(mABsq,mACsq);
       if(resnum<2) hel = CosHelD0pi2(mABsq,mBCsq,pq);// D pi resonance
  else if(resnum<5) hel = CosHelpipi2(mBCsq,mABsq,pq);// pi pi resonance
  switch(resnum){
  case 0: return DstarAmp(  mABsq,hel,pq); // D*
  case 1: return DdstarAmp( mABsq,hel,pq); // D2
  case 2: return DzeroAmp(  mABsq       ); // D0
  case 3: return FullRhoAmp(mBCsq,hel,pq); // rho + omega + rho'
  case 4: return f2Amp(     mBCsq,hel,pq); // f2
  case 5: return f0Amp(     mBCsq       ); // f0
  case 6: return f980Amp(   mBCsq       ); // f980
  case 7: return f1370Amp(  mBCsq       ); // f1370
  default: cout << "KAmp: wrong resnum " << resnum << endl;
  }
  return 0;
}

void KAmp::GetResVals(vectcd& resv, cdouble& mABsq,cdouble& mACsq) const{
  resv.clear();
  cdouble mBCsq = GetmBCsq(mABsq,mACsq);
  double pq2;
  double pq3;
  cdouble h2 = CosHelD0pi2(mABsq,mBCsq,pq2);
  cdouble h3 = CosHelpipi2(mBCsq,mACsq,pq3);
  resv.push_back(DstarAmp( mABsq,h2,pq2));// D* amplitude
  resv.push_back(DdstarAmp(mABsq,h2,pq2));// D** amplitude
  resv.push_back(DzeroAmp( mABsq));       // D0* amplitude
  resv.push_back(rhoAmp(   mBCsq,h3,pq3));// rho(770) amplitude
  resv.push_back(omegaAmp( mBCsq,h3,pq3));// omega amplitude
  resv.push_back(rhopAmp(  mBCsq,h3,pq3));// rho' amplitude
  resv.push_back(f2Amp(    mBCsq,h3,pq3));// f2 amplitude
  resv.push_back(f0Amp(    mBCsq));       // f0 amplitude
  resv.push_back(f980Amp(  mBCsq));       // f(980) amplitude
  resv.push_back(f1370Amp( mBCsq));       // f(1370) amplitude
}

void KAmp::init(void){
  SetCoefficients({m_dm1,m_dm2,m_dm0,m_dmrh,m_dmf2,m_dmf0,m_dmbs,m_d3});
  SetAmpSignature({1,1,1,3,1,1,1,1});

  const double nsig = 0.5;
  vectd ledge = {
      pow(2.01-0.05,2),
      pow(m_ams2-nsig*m_gms2,2),
      pow(m_ams0-nsig*m_gms0,2),
      pow(m_amrh,2),
      pow(m_amf2-nsig*m_gmf2,2),
      pow(m_amf0-m_gmf0,2),
      pow(m_amf098-nsig*m_gmf098,2),
      pow(m_am3-nsig*m_gm3,2)
  };
  vectd redge = {
      pow(2.01+nsig*0.05,2),
      pow(m_ams2+nsig*m_gms2,2),
      pow(m_ams0+nsig*m_gms0,2),
      pow(m_amrh1,2),
      pow(m_amf2+nsig*m_gmf2,2),
      pow(m_amf0+m_gmf0,2),
      pow(m_amf098+nsig*m_gmf098,2),
      pow(m_am3+nsig*m_gm3,2)
  };
  vecti types {1,1,1,3,3,3,3,3};
  SetResAreas(ledge,redge,types);
}

double KAmp::Norm(void) const{ return NormWithCache();}

void KAmp::SetParams(const vectd& pars){
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

//double KAmp::ang1(cdouble& h,cdouble& pq){ return 16.*pq*h;}
//double KAmp::ang2(cdouble& h,cdouble& pq){ return 4.*pow(pq,2)*(pow(h,2)-1./3.);}
double KAmp::ang1(cdouble& h,cdouble& pq){ return pq*h;}
double KAmp::ang2(cdouble& h,cdouble& pq){ return pow(pq,2)*(pow(h,2)-1./3.);}

//compld KAmp::btodpp(cdouble& qp,cdouble& qn) const {
//  vectcd coefv;
//  GetResAmp(coefv,qp,qn);
//  compld res = 0;
//  for(unsigned i=0; i<coefv.size(); i++) res += m_ampl[i]*coefv[i];
//  return res;
//}

void KAmp::PrintSummary(void){
  cout << "*** The Kuzmin amplitude ***" << endl;
  cout << "Masses: B0 " << mM() << ", D0 " << mA() << ", pi+ " << mB() << endl;
  cout << "D*      amp: " << m_dm1  << endl;
  cout << "D2      amp: " << m_dm2  << endl;
  cout << "D0      amp: " << m_dm0  << endl;
  cout << "rho     amp: " << m_dmrh << ", anrh: " << anrh(m_amrh,m_gmrh) << endl;
  cout << "omega   amp: " << m_dmrh << endl;
  cout << "rho'    amp: " << m_dmrh << endl;
  cout << "f2      amp: " << m_dmf2 << ", anf2: " << anf2(m_amf2,m_gmf2) << endl;
  cout << "f(980)  amp: " << m_dmbs << endl;
  cout << "f0(550) amp: " << m_dmf0 << ", anf0: " << anf0(m_amf0,m_gmf0) << endl;
  cout << "f(1370) amp: " << m_d3   << ", anf0: " << anf0(m_am3, m_gm3)  << endl;
  return;
}

//compld KAmp::Amp(cdouble& mAB, cdouble& mAC) const{return btodpp(mAB,mAC);}

double KAmp::GetmBCsq(cdouble& mABsq,cdouble& mACsq){return m_fm2_2pm2_md2-mABsq-mACsq;}

// M -> (R -> AB) C. Returns A energy in the R frame
double KAmp::eA(cdouble& mAsq, cdouble& mBsq, cdouble& mABsq){
  return (mAsq + mABsq - mBsq)/(2.*sqrt(mABsq));
}
// M -> (R -> AB) C. Returns A energy in the R frame
double KAmp::eB(cdouble& mAsq, cdouble& mBsq, cdouble& mABsq){
  return (mBsq + mABsq - mAsq)/(2.*sqrt(mABsq));
}
// M -> (R -> AB) C. Returns C energy in the R frame
double KAmp::eC(cdouble& mMsq, cdouble& mCsq, cdouble& mABsq){
  return (mMsq - mABsq - mCsq)/(2.*sqrt(mABsq));
}
// M -> (R -> AB) C. Returns C energy in the M frame
double KAmp::eC0(cdouble& mMsq, cdouble& mCsq, cdouble& mABsq){
  return (mMsq + mCsq - mABsq)/(2.*sqrt(mMsq));
}

// A = D0, B = pi, C = pi
double KAmp::CosHelD0pi(cdouble& mDpipsq, cdouble& mDpimsq, double& pq){
  return MyCosHelAB(m_fm2,m_dstm2,m_pm2,m_pm2,mDpipsq,mDpimsq,pq);
}
// A = pi, B = pi, C = D0
double KAmp::CosHelpipi(cdouble& mpipisq, cdouble& mDpipsq, double& pq){
  return MyCosHelAB(m_fm2,m_pm2,m_pm2,m_dstm2,mpipisq,mDpipsq,pq);
}
// A = D0, B = pi, C = pi
double KAmp::CosHelD0pi2(cdouble& mDpipsq, cdouble& mpipimsq, double& pq){
  return MyCosHelAB2(m_fm2,m_dstm2,m_pm2,m_pm2,mDpipsq,mpipimsq,pq);
}
// A = pi, B = pi, C = D0
double KAmp::CosHelpipi2(cdouble& mpipisq, cdouble& mDpimsq, double& pq){
  return MyCosHelAB2(m_fm2,m_pm2,m_pm2,m_dstm2,mpipisq,mDpimsq,pq);
}

double KAmp::MyCosHelAB(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq, cdouble& mACsq, double& pq){
  cdouble enA   = eA(mAsq,mBsq,mABsq);
  cdouble moAsq = pow(enA,2) - mAsq;
  cdouble moA   = moAsq>0 ? sqrt(moAsq) : 0;
  cdouble enC   = eC(mMsq,mCsq,mABsq);
  cdouble moCsq = pow(enC,2) - mCsq;
  cdouble moC   = moCsq > 0 ? sqrt(moCsq) : 0;
  pq = moA*moC;
  return (mAsq + mCsq + 2.*enA*enC - mACsq)/(2.*moA*moC);
}

double KAmp::MyCosHelAB2(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq, cdouble& mBCsq, double& pq){
  cdouble enB   = eB(mAsq,mBsq,mABsq);
  cdouble moBsq = pow(enB,2) - mBsq;
  if(moBsq<0){ pq = 0; return 0;}
//    cout << "MyCosHelAB2: moBsq " << moBsq << " " << enB << " " << mBsq << " " << mABsq << " " << mAsq << endl;
  cdouble moB   = sqrt(moBsq);

  cdouble enC   = eC(mMsq,mCsq,mABsq);
  cdouble moCsq = pow(enC,2) - mCsq;
  if(moCsq<0){ pq = 0; return 0;}
//    cout << "MyCosHelAB2: moCsq " << moCsq << " " << enC << " " << mCsq << endl;
  cdouble moC   = sqrt(moCsq);

//  cdouble enC0  = eC0(mMsq,mCsq,mABsq);
//  cdouble moC0sq = pow(enC0,2) - mCsq;
//  if(moC0sq<0){ pq = 0; return 0;}
//  cdouble moC0 = sqrt(moC0sq);

//  pq = moB*moC0;
  pq = moB*moC;
  return (mBsq + mCsq + 2.*enB*enC - mBCsq)/(2.*moB*moC);
}

double KAmp::CosHelAB(cdouble& mMo,cdouble& mA,cdouble& mB,cdouble& mC,cdouble& mABsq){
  cdouble mABsqmin = pow(mA+mB,2);
  cdouble mABsqmax = pow(mMo-mC,2);
  return (mABsqmax+mABsqmin-2.*mABsq)/(mABsqmax-mABsqmin);
}
// M -> (R -> AB) C. Returns Dot of p(A) and p(C) 3-vectors in the R frame
double KAmp::DotAC(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq, cdouble& mACsq){
  cdouble enA = eA(mAsq,mBsq,mABsq);
  cdouble enC = eC(mMsq,mBsq,mABsq);
  return 0.5*(mAsq + mCsq + 2.*enA*enC - mACsq);
}

compld KAmp::amrs0(cdouble& q2,cdouble& am,cdouble& gm){
  cdouble am2  = pow(am,2);
  cdouble ar   = q2-am2;
  cdouble q    = sqrt(q2);
  // pion energy
  cdouble epi0 = (am2+m_pm2-m_dstm2)/(2.*am);
  cdouble epi  = (q2 +m_pm2-m_dstm2)/(2.*q);
  // pion momentum
  cdouble ppi0 = sqrt(pow(epi0,2)-m_pm2);
  cdouble ppi  = sqrt(pow(epi ,2)-m_pm2);
  // resonance energy
  cdouble eres0 = (m_fm2+am2-m_pm2)/(2*m_fm);
  if(pow(eres0,2)<am2){
    cout << "amrs0: eres0<am2 " << eres0 << " " << am << endl;
    return 0;
  }
  cdouble gr = (ppi/ppi0)*(gm*am2/q);
  return 1./compld(ar,gr);
}
// D0 pi resonances
compld KAmp::amrs1(cdouble& q2){
  // D* veto
  if(fabs(q2-4.04)<0.01) return compld(0,0);
  cdouble am    = 2.01;
  cdouble gm    = 0.0001;
  cdouble am2   = pow(am,2);
  cdouble q     = sqrt(q2);
  // pion energy in the resonance frame calculated with nominal resonance mass
  cdouble epi0  = (am2+m_pm2-m_dstm2)/(2.*am);
  cdouble epi   = (q2 +m_pm2-m_dstm2)/(2.*q);
  // pion momentum in the resonance frame calculated with nominal resonance mass
  cdouble ppi0  = sqrt(pow(epi0,2)-m_pm2);
  cdouble ppi   = sqrt(pow(epi ,2)-m_pm2);
  // resonance energy in the B frame calculated with nominal resonance mass
  cdouble eres0 = (m_fm2+am2-m_pm2)/(2.*m_fm);
  cdouble eres  = (m_fm2+q2 -m_pm2)/(2.*m_fm);
  // resonance momentum
  cdouble pres0 = sqrt(pow(eres0,2)-am2);
  cdouble pres  = sqrt(pow(eres ,2)-am2);
  // B0 formfactor
  cdouble rr_b  = pres *m_ffr;
  cdouble r0_b  = pres0*m_ffr;
  cdouble fg    = sqrt((1+pow(r0_b,2))/(1+pow(rr_b,2)));
  // resonance formfactor
  cdouble rr_r  = ppi *m_ffdel;
  cdouble r0_r  = ppi0*m_ffdel;
  cdouble f     = exp(-rr_r+r0_r);
  // full Breit-Wigner
  cdouble ar = q2-am2;
  cdouble gr = (gm*am2/q)*pow(ppi/ppi0,3)*pow(f,2);
  return f*fg/(compld(ar,gr)*3.53438);
}

compld KAmp::amrs2(cdouble& q2,cdouble& am,cdouble& gm){
  cdouble am2 = pow(am,2);
  cdouble q   = sqrt(q2);
  // pion energy
  cdouble epi0  = (am2+m_pm2-m_dstm2)/(2.*am);
  cdouble epi   = (q2 +m_pm2-m_dstm2)/(2.*q);
  // pion momentum
  cdouble ppi0  = sqrt(pow(epi0,2)-m_pm2);
  cdouble ppi   = sqrt(pow(epi ,2)-m_pm2);
  // resonance energy
  cdouble eres0 = (m_fm2+am2-m_pm2)/(2.*m_fm);
  cdouble eres  = (m_fm2+q2-m_pm2)/(2.*m_fm);
  // resonance momentum
  cdouble pres0 = sqrt(pow(eres0,2)-am2);
  cdouble pres  = sqrt(pow(eres ,2)-am2);
  // B0 formfactor
  cdouble rr_b  = pres *m_ffr;
  cdouble r0_b  = pres0*m_ffr;
  cdouble fg    = sqrt((9+3*pow(r0_b,2)+pow(r0_b,4))/(9+3*pow(rr_b,2)+pow(rr_b,4)));
  // resonance formfactor
  cdouble rr_r  = ppi*m_ffr;
  cdouble r0_r  = ppi0 *m_ffr;
  cdouble ff    = sqrt((9.+3.*pow(r0_r,2)+pow(r0_r,4))/(9.+3.*pow(rr_r,2)+pow(rr_r,4)));
  // full Breit-Wigner
  cdouble ar = q2-am2;
  cdouble gr = (gm*am2/q)*pow(ppi/ppi0,5)*pow(ff,2);
  return ff*fg/compld(ar,gr);
}

// pi+ pi- resonances
compld KAmp::amrs0h(cdouble& q2,cdouble& am,cdouble& gm){
  cdouble am2  = pow(am,2);
  cdouble q    = sqrt(q2);
  // pion energy
  cdouble epi0 = 0.5*am;
  cdouble epi  = 0.5*q;
  // pion momentum
  cdouble ppi0 = sqrt(pow(epi0,2)-m_pm2);
  cdouble ppi  = sqrt(pow(epi ,2)-m_pm2);
  // resonance energy
  cdouble eres0 = (m_fm2+am2-m_dstm2)/(2.*m_fm);
  if(pow(eres0,2)<am2){
     cout << "amrs0h: eres0<am2 " << eres0 << " " << am << endl;
     return 0;
  }
  cdouble ar = q2-am2;
  cdouble gr = (gm*am2/q)*(ppi/ppi0);
  return 1./compld(ar,gr);
}

compld KAmp::amrs1h(cdouble& q2,cdouble& am,cdouble& gm){
  cdouble am2  = pow(am,2);
  cdouble q    = sqrt(q2);
  // pion energy
  cdouble epi0 = 0.5*am;
  cdouble epi  = 0.5*q;
  // pion momentum
  cdouble ppi0 = sqrt(pow(epi0,2)-m_pm2);
  cdouble ppi  = sqrt(pow(epi ,2)-m_pm2);
  // resonance energy
  cdouble eres0 = (m_fm2+am2-m_dstm2)/(2*m_fm);
  if(pow(eres0,2)<am2){
     cout << "amrs1h: eres0<am2 " << eres0 << " " << am << endl;
     return 0;
  }
  double hwd0,hwd;
  cdouble ar  = q2-am2;
  cdouble gr  = (gm*am2/q)*pow(ppi/ppi0,3);
  cdouble hw  = hwrho(q2,hwd);
  cdouble hw0 = hwrho(am2,hwd0);
  cdouble dm  = ((hw-hw0)*pow(ppi,2)-ar*hwd0*pow(ppi0,2))*am2*gm/pow(ppi0,3);
  return 1./compld(ar-dm,gr);
}

double KAmp::hwrho(cdouble& s, double& hwd){
  cdouble y = sqrt(1.-0.07795/s);
  cdouble w = log((1.+y)/(1-y));
  hwd = (0.5*w*(1.-pow(y,2))/y+1.)/s*M_PI_2;
  return w*y*M_PI_2;
}

compld KAmp::amrs2h(cdouble& q2,cdouble& am,cdouble& gm){
  cdouble am2  = pow(am,2);
  cdouble q    = sqrt(q2);
  // pion energy
  cdouble epi0 = 0.5*am;
  cdouble epi  = 0.5*q;
  // pion momentum
  cdouble ppi0 = sqrt(pow(epi0,2)-m_pm2);
  cdouble ppi  = sqrt(pow(epi ,2)-m_pm2);
  // resonance energy
  cdouble eres0 = 0.5*(m_fm2+am2-m_dstm2)/m_fm;
  if(pow(eres0,2)<am2){
    cout << "amrs1h: eres0<am " << eres0 << " " << am << endl;
    return 0;
  }
  cdouble ar = q2-am2;
  cdouble gr = (gm*am2/q)*pow(ppi/ppi0,5);
  return 1./compld(ar,gr);
}

// omega resonanse
compld KAmp::amrs1hm(cdouble& q2){
  cdouble am  = 0.78257;
  cdouble am2 = pow(am,2);
  // pion energy
  cdouble epi0  = 0.5*(m_fm2+am2-m_dstm2)/m_fm;
  if(pow(epi0,2)<am2){
    cout << "amrs1hm: epi0<am2 " << epi0 << " " << am << endl;
    return 0;
  }
  cdouble ar = q2-am2;
  cdouble gr = sqrt(q2)*0.00849;
  return 1./compld(ar,gr);
}

// amplitude wrappers
double KAmp::anf0(cdouble& am, cdouble& ag){
  const vectd pars = {  -69.398,
                        608.83,
                      -1152.1,
                        940.13,
                       -198.62,
                       -148.03,
                         66.298};
  cdouble poly = pars[0]+pars[1]*am+
                      pars[2]*pow(am,2)+
                      pars[3]*pow(am,3)+
                      pars[4]*pow(am,4)+
                      pars[5]*pow(am,5)+
                      pars[6]*pow(am,6);
  cdouble result = poly/(pow(am,0.6)*pow(ag,1.25));
  if(isnan(result)){
    cout << "anf0: bad result" << endl;
    return 1;
  }
  return result;
}

double KAmp::anf2(cdouble& am, cdouble& ag){
  const vectd pars =  {4635.9,
                        -21.723,
                        388.72,
                        -94.873,
                       1162.3,
                      -4993.5};
  cdouble poly = pars[0]+pars[1]*am+pars[2]*ag+
                      pars[3]*pow(am,2)+pars[4]*am*ag+pars[5]*pow(ag,2);
  return poly/pow(am,1.55)/pow(ag,1.095);
}

double KAmp::anrh(cdouble& am, cdouble& ag){
  return 546.8/pow(am,0.6)/pow(ag,1.12);
}

double KAmp::an2( cdouble& am, cdouble& ag){
  const vectd pars = { 0.18907E-02,
                      -0.50839E-05,
                       0.16955E-02,
                      -0.19040E-05,
                       0.69806E-03,
                      -0.59716E-01};
  cdouble poly = pars[0]+pars[1]*am+pars[2]*ag+
                      pars[3]*pow(am,2)+pars[4]*am*ag+pars[5]*pow(ag,2);
  return poly*pow(am,10.5)/pow(ag,1.02);
}

double KAmp::an0( cdouble& am, cdouble& ag){
  const vectd pars = {-9.5637,
                       8.7421,
                      -0.28803,
                      -1.8988,
                       0.40894,
                      -1.5582,
                      -0.84441E-03,
                       0.29120E-01,
                      -0.32185,
                       1.9611};
  cdouble poly = pars[0]+pars[1]*am+pars[2]*ag+
                      pars[3]*pow(am,2)+pars[4]*am*ag+pars[5]*pow(ag,2)+
                      pars[6]*pow(am,3)+pars[7]*pow(am,2)*ag+pars[8]*am*pow(ag,2)+pars[9]*pow(ag,3);
  return poly*pow(am,3)/pow(ag,1.15);
}

// Row amplitudes and phases
cdouble KAmp::m_a0  = 0.035616; // par 1      D0 'raw' amplitude
cdouble KAmp::m_f0  =-3.86950;  // par 3      D0 phase
cdouble KAmp::m_a1  = 0.077521; // par 2      D* 'raw' amplitude
cdouble KAmp::m_f1  = 2.75360;  // par 4      D* phase
cdouble KAmp::m_ab  = 0;        // par 10     f(980) 'raw' amplitude
cdouble KAmp::m_fb  = 1.0030;   // par 11     f(980) phase
cdouble KAmp::m_a3  = 1.961500; // par 24     f(1370) 'raw' amplitude
cdouble KAmp::m_f3  = 0.81462;  // par 25     f(1370) phase
cdouble KAmp::m_ar  = 0.428140; // par 14
cdouble KAmp::m_fr  = 1.43780;  // par 15
cdouble KAmp::m_af0 = 0.054742; // par 18
cdouble KAmp::m_ff0 =-1.72320;  // par 19
cdouble KAmp::m_af2 = 0.081368; // par 22
cdouble KAmp::m_ff2 = 2.71440;  // par 23

// Masses and widths
cdouble KAmp::m_fm     = 5.279;
cdouble KAmp::m_pm     = 0.139568;
cdouble KAmp::m_dstm   = 1.865;
cdouble KAmp::m_ams0   = 2.3080;   // par 7
cdouble KAmp::m_gms0   = 0.27611;  // par 8
cdouble KAmp::m_ams2   = 2.4677;   // par 5
cdouble KAmp::m_gms2   = 0.056002; // par 6
cdouble KAmp::m_amrh   = 0.77560;  // par 12
cdouble KAmp::m_gmrh   = 0.15000;  // par 13
cdouble KAmp::m_amrh1  = 1.465;
cdouble KAmp::m_gmrh1  = 0.31;
cdouble KAmp::m_amf0   = 0.51300;  // par 16
cdouble KAmp::m_gmf0   = 0.33500;  // par 17
cdouble KAmp::m_amf098 = 0.978;    // fixed
cdouble KAmp::m_gmf098 = 0.0440;   // fixed
cdouble KAmp::m_am3    = 1.434;
cdouble KAmp::m_gm3    = 0.173;
cdouble KAmp::m_amf2   = 1.2750;   // par 20
cdouble KAmp::m_gmf2   = 0.18500;  // par 21
cdouble KAmp::m_amom   = 0.0272;   // par 30

cdouble KAmp::m_fm2_2pm2_md2 = m_fm2+2*m_pm2+m_dstm2;

// Auxiliary parameters
cdouble KAmp::m_fm2    = pow(KAmp::m_fm,2);
cdouble KAmp::m_pm2    = pow(KAmp::m_pm,2);
cdouble KAmp::m_dstm2  = pow(KAmp::m_dstm,2);
cdouble KAmp::m_ddd    = get_ddd();

// Complex amplitudes
const compld KAmp::m_dm2  = get_dm2();
const compld KAmp::m_dm0  = get_dm0();
const compld KAmp::m_dm1  = get_dm1();
const compld KAmp::m_dmrh = get_dmrh();
const compld KAmp::m_dmf2 = get_dmf2();
const compld KAmp::m_dmf0 = get_dmf0();
const compld KAmp::m_dmbs = get_dmbs();
const compld KAmp::m_d3   = get_d3();

// Amplitude wrappers
double KAmp::get_ddd(void)  { return 1. - m_a0 + m_a1 + m_ab + m_a3 + m_ar + m_af0 + m_af2;}
compld KAmp::get_dm0(void)  { return exp(      m_f0  *ione)*sqrt(m_a0/an0(m_ams0,m_gms0));}
compld KAmp::get_dm1(void)  { return exp(      m_f1  *ione)*sqrt(m_a1);}
double KAmp::get_dm2(void)  { return                        sqrt(get_ddd()/an2(m_ams2,m_gms2));}
compld KAmp::get_dmbs(void) { return exp(      m_fb  *ione)*sqrt(m_ab/anf0(m_amf098,m_gmf098));}
compld KAmp::get_d3(void)   { return exp(      m_f3  *ione)*sqrt(m_a3/anf0(m_am3,m_gm3));}
compld KAmp::get_dmrh(void) { return exp(      m_fr  *ione)*sqrt(m_ar/anrh(m_amrh,m_gmrh));}
compld KAmp::get_dmf0(void) { return exp((m_fr+m_ff0)*ione)*sqrt(m_af0/anf0(m_amf0,m_gmf0));}
compld KAmp::get_dmf2(void) { return exp((m_fr+m_ff2)*ione)*sqrt(m_af2/anf2(m_amf2,m_gmf2));}