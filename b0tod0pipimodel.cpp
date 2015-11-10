#include "b0tod0pipimodel.h"

B0toD0pipiModel::B0toD0pipiModel(void):
B0toD0pipiModel(5.27958,1.86484,0.13957018)
{
}

B0toD0pipiModel::B0toD0pipiModel(const double& mB, const double& mD, const double& mpi):
SymDalitzModel(mB,mD,mpi,-M_PI/8.,15.*M_PI/8.)
{
// ** A. Kuzmin et al. (Belle Collaboration) Phys. Rev. D 76, 012006 – Published 30 July 2007 **
  const double radtodeg = EvtConst::radToDegrees;
  EvtVector4R p4_p,moms1,moms2,moms3;
  AddRes(new EvtResonance2(p4_p,moms1,moms2, 2.15, 0.00*radtodeg, 0.0496,  2.4657,  2));//D2*
  AddRes(new EvtResonance2(p4_p,moms1,moms2, 0.60,-3.00*radtodeg, 0.2760,  2.3080,  0));//D0*
  AddRes(new EvtResonance2(p4_p,moms1,moms2, 0.88,-2.62*radtodeg, 83.4e-6, 2.01027, 1));//Dv*

  AddRes(new EvtResonance2(p4_p,moms2,moms3, 3.19, 2.25*radtodeg, 0.144,   0.7756,  1));//rho
  AddRes(new EvtResonance2(p4_p,moms2,moms3, 0.68, 2.97*radtodeg, 0.185,   1.275,   2));//f2
  AddRes(new EvtResonance2(p4_p,moms2,moms3, 0.68,-0.44*radtodeg, 0.335,   0.513,   0));//f0(600)
  AddRes(new EvtResonance2(p4_p,moms2,moms3, 0.08,-2.48*radtodeg, 0.044,   0.978,   0));//f0(980)
  AddRes(new EvtResonance2(p4_p,moms2,moms3, 0.21,-1.52*radtodeg, 0.173,   1.434,   0));//f0(1470)
}

EvtComplex B0toD0pipiModel::Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3){
  return amp_BelleKuzmin(p4_p,moms1,moms2,moms3);
}

EvtComplex B0toD0pipiModel::amp_BelleKuzmin(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3){
  EvtComplex amp(0.0);
  SetMomenta(0,p4_p,moms1,moms2);//D2*
  SetMomenta(1,p4_p,moms1,moms2);//D0*
  SetMomenta(2,p4_p,moms1,moms2);//Dv*

  SetMomenta(3,p4_p,moms2,moms3);//rho
  SetMomenta(4,p4_p,moms2,moms3);//f2
  SetMomenta(5,p4_p,moms2,moms3);//f0(600)
  SetMomenta(6,p4_p,moms2,moms3);//f0(980)
  SetMomenta(7,p4_p,moms2,moms3);//f0(1470)
  for(int i=0; i<ResNum(); i++){ amp += const_cast<EvtResonance2*>(Res(i))->resAmpl();}
  if(std::isnan(real(amp)) || std::isnan(imag(amp))){
    cout << "amp_BelleKuzmin: (" << real(amp) << "," << imag(amp) << ")" << endl;
    cout << " p4_p  " << p4_p << endl;
    cout << " moms1 " << moms1 << endl;
    cout << " moms2 " << moms2 << endl;
    cout << " moms3 " << moms3 << endl;
  }
  return amp;
}
