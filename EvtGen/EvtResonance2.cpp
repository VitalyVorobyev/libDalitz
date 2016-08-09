#include "EvtResonance2.h"

#include <cmath>
#include "EvtVector4R.h"
#include "EvtKine.h"
#include "EvtComplex.h"
//#include "EvtReport.h"
#include "EvtConst.h"

EvtResonance2::~EvtResonance2(){}

//operator
EvtResonance2& EvtResonance2::operator = ( const EvtResonance2  &n)
{
  if ( &n == this ) return *this;
  _p4_p = n._p4_p;
  _p4_d1 = n._p4_d1;
  _p4_d2 = n._p4_d2;
  _ampl = n._ampl;
  _theta = n._theta;
  _gamma = n._gamma;
  _spin = n._spin;
  _bwm = n._bwm;
   return  *this;
}

//constructor

EvtResonance2::EvtResonance2(const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
                           const  EvtVector4R& p4_d2, double ampl,
                           double theta, double gamma, double bwm, int spin):
  _p4_p(p4_p),_p4_d1(p4_d1), _p4_d2(p4_d2),_ampl(ampl), _theta(theta/EvtConst::radToDegrees),
  _gamma(gamma), _bwm(bwm), _spin(spin) {}

//amplitude function
EvtComplex EvtResonance2::resAmpl(){
//  double pi180inv = 1.0/EvtConst::radToDegrees;
  EvtComplex ampl;
  EvtVector4R  p4_d3 = _p4_p-_p4_d1-_p4_d2;

  //get cos of the angle between the daughters from their 4-momenta
  //and the 4-momentum of the parent

  //in general, EvtDecayAngle(parent, part1+part2, part1) gives the angle
  //the missing particle (not listed in the arguments) makes
  //with part2 in the rest frame of both
  //listed particles (12)

  //angle 3 makes with 2 in rest frame of 12 (CS3)
//  double cos_phi_0 = EvtDecayAngle(_p4_p, _p4_d1+_p4_d2, _p4_d1);
  //angle 3 makes with 1 in 12 is, of course, -cos_phi_0

  //first compute several quantities...follow CLEO preprint 00-23

  double mAB=(_p4_d1+_p4_d2).mass();
  double mBC=(_p4_d2+p4_d3).mass();
  double mAC=(_p4_d1+p4_d3).mass();
  double mA=_p4_d1.mass();
  double mB=_p4_d2.mass();
  double mD=_p4_p.mass();
  double mC=p4_d3.mass();

  double mR=_bwm;
  double gammaR=_gamma;
  double pAB=sqrt( (((mAB*mAB-mA*mA-mB*mB)*(mAB*mAB-mA*mA-mB*mB)/4.0) -
                    mA*mA*mB*mB)/(mAB*mAB));
  double pR=sqrt( (((mR*mR-mA*mA-mB*mB)*(mR*mR-mA*mA-mB*mB)/4.0) -
                   mA*mA*mB*mB)/(mR*mR));

  double pD= (((mD*mD-mR*mR-mC*mC)*(mD*mD-mR*mR-mC*mC)/4.0) -
                   mR*mR*mC*mC)/(mD*mD);
  if ( pD>0 ) { pD=sqrt(pD); } else {pD=0;}
  double pDAB=sqrt( (((mD*mD-mAB*mAB-mC*mC)*(mD*mD-mAB*mAB-mC*mC)/4.0) -
                   mAB*mAB*mC*mC)/(mD*mD));

  //    report(INFO,"EvtGen") << mAB<<" "<< mBC<<" "<< mAC<<" "<< mA<<" "<< mB<<" "<< mC<<" "
  //     << mD<<" "<< mR<<" "<< gammaR<<" "<< pAB<<" "<< pR<<" "<< pD<<" "<<pDAB<<std::endl;

  double fR=1;
  double fD=1;
  int power;
  switch (_spin) {
  case 0:
    fR=1.0;
    fD=1.0;
    power=1;
    //report(INFO,"EvtGen") << "fR="<<fR<<" fD="<<fD<<std::endl;
    break;
  case 1:
    fR=sqrt(1.0+1.5*1.5*pR*pR)/sqrt(1.0+1.5*1.5*pAB*pAB);
    fD=sqrt(1.0+5.0*5.0*pD*pD)/sqrt(1.0+5.0*5.0*pDAB*pDAB);
    //report(INFO,"EvtGen") << "fR="<<fR<<" fD="<<fD<<std::endl;
    power=3;
    break;
  case 2:
    fR = sqrt( (9+3*pow((1.5*pR),2)+pow((1.5*pR),4))/(9+3*pow((1.5*pAB),2)+pow((1.5*pAB),4)) );
    fD = sqrt( (9+3*pow((5.0*pD),2)+pow((5.0*pD),4))/(9+3*pow((5.0*pDAB),2)+pow((5.0*pDAB),4)) );
    power=5;
    //report(INFO,"EvtGen") << "fR="<<fR<<" fD="<<fD<<std::endl;
    break;
  default:
//    report(INFO,"EvtGen") << "Incorrect spin in EvtResonance22.cc\n";
    std::cout << "Incorrect spin in EvtResonance2.cpp\n";
  }

  double gammaAB= gammaR*pow(pAB/pR,power)*(mR/mAB)*fR*fR;
  //report(INFO,"EvtGen") << gammaAB<<std::endl;
  switch (_spin) {
  case 0:
//    ampl=_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
    ampl=_ampl*EvtComplex(cos(_theta),sin(_theta))*
          fR*fD/(mR*mR-mAB*mAB-EvtComplex(0.0,mR*gammaAB));
    break;
  case 1:
//    ampl=_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
    ampl=_ampl*EvtComplex(cos(_theta),sin(_theta))*
      (fR*fD*(mAC*mAC-mBC*mBC+((mD*mD-mC*mC)*(mB*mB-mA*mA)/(mR*mR)))/
       (mR*mR-mAB*mAB-EvtComplex(0.0,mR*gammaAB)));
    break;
  case 2:
  //     ampl=_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
  //      fR*fD/(mR*mR-mAB*mAB-EvtComplex(0.0,mR*gammaAB))*
  //        (pow((mBC*mBC-mAC*mAC+(mD*mD-mC*mC)*(mA*mA-mB*mB)/(mAB*mAB)),2)-
  //      (1.0/3.0)*(mAB*mAB-2*mD*mD-2*mC*mC+pow((mD*mD- mC*mC)/mAB, 2))*
  //      (mAB*mAB-2*mA*mA-2*mB*mB+pow((mA*mA-mB*mB)/mAB,2)));
//    ampl=_ampl*EvtComplex(cos(_theta*pi180inv),sin(_theta*pi180inv))*
    ampl=_ampl*EvtComplex(cos(_theta),sin(_theta))*
         fR*fD/(mR*mR-mAB*mAB-EvtComplex(0.0,mR*gammaAB))*
         (pow((mBC*mBC-mAC*mAC+(mD*mD-mC*mC)*(mA*mA-mB*mB)/(mR*mR)),2)-
         (1.0/3.0)*(mAB*mAB-2*mD*mD-2*mC*mC+pow((mD*mD- mC*mC)/mR, 2))*
         (mAB*mAB-2*mA*mA-2*mB*mB+pow((mA*mA-mB*mB)/mR,2)));
    break;
  default:
//    report(INFO,"EvtGen") << "Incorrect spin in EvtResonance22.cc\n";
    std::cout << "Incorrect spin in EvtResonance2.cpp\n";
  }

  //report(INFO,"EvtGen") <<"The amplitude is "<<ampl<<std::endl;
  return ampl;
}
