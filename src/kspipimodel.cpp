#include "kspipimodel.h"

KspipiModel::KspipiModel(void):
KspipiModel(1.86484,0.497614,0.13957018)
{
}
KspipiModel::KspipiModel(const double &md, const double &mks, const double& mpi):
SymDalitzModel(md,mks,mpi,-M_PI/8.,15.*M_PI/8.)
{
  SetABaxis("m_{+}^{2}, GeV^{2}/c^{4}");
  SetACaxis("m_{-}^{2}, GeV^{2}/c^{4}");
  SetBCaxis("m_{#pi#pi}^{2}, GeV^{2}/c^{4}");
  // ** A. Poluektov et al. Phys. Rev. D 81, 112002 – Published 16 June 2010 **
  EvtVector4R p4_p,moms1,moms2,moms3;
  AddRes(new EvtResonance2(p4_p,moms1,moms2, 1.638, 133.2, 0.0484, 0.8937, 1));//K*(892)
  AddRes(new EvtResonance2(p4_p,moms1,moms2, 2.210, 358.9, 0.294 , 1.412 , 0));//K0*(1430)
  AddRes(new EvtResonance2(p4_p,moms1,moms2, 0.890, 314.8, 0.0985, 1.4256, 2));//K2*(1430)
  AddRes(new EvtResonance2(p4_p,moms1,moms2, 0.880,  82.0, 0.322 , 1.717 , 1));//K*(1680)
  AddRes(new EvtResonance2(p4_p,moms1,moms2, 0.650, 120.0, 0.232 , 1.414 , 1));//K*(1410)

  AddRes(new EvtResonance2(p4_p,moms1,moms3, 0.149, 325.4, 0.0508, .89166, 1));//DCS K*(892)
  AddRes(new EvtResonance2(p4_p,moms1,moms3, 0.360,  87.0, 0.294 , 1.412 , 0));//DCS K0*(1430)
  AddRes(new EvtResonance2(p4_p,moms1,moms3, 0.230, 275.0, 0.0985, 1.4256, 2));//DCS K2*(1430)
  AddRes(new EvtResonance2(p4_p,moms1,moms3, 2.100, 130.0, 0.322 , 1.717 , 1));//DCS K*(1680)
  AddRes(new EvtResonance2(p4_p,moms1,moms3, 0.420, 253.0, 0.232 , 1.414 , 1));//DCS K*(1410)

  AddRes(new EvtResonance2(p4_p,moms3,moms2, 1.000,   0.0, 0.1490, 0.7717, 1));//Rho
  AddRes(new EvtResonance2(p4_p,moms3,moms2, .0343, 112.0, .00849, .78265, 1));//Omega
  AddRes(new EvtResonance2(p4_p,moms3,moms2, 0.385, 207.3, 0.05  , 0.977,  0));//f0(980)
  AddRes(new EvtResonance2(p4_p,moms3,moms2, 1.250,  69.0, 0.272 , 1.31 ,  0));//f0(1370)
  AddRes(new EvtResonance2(p4_p,moms3,moms2, 1.440, 342.9, 0.1851, 1.2754, 2));//f2(1270)
  AddRes(new EvtResonance2(p4_p,moms3,moms2, 0.490,  64.0, 0.400 , 1.465,  1));//Rho(1450)
  AddRes(new EvtResonance2(p4_p,moms3,moms2, 1.560, 214.0, 0.453 , 0.522,  0));//sigma1
  AddRes(new EvtResonance2(p4_p,moms3,moms2, 0.200, 212.0, 0.088 , 1.033,  0));//sigma2
}

EvtComplex KspipiModel::Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3){
  return amp_Belle2010(p4_p,moms1,moms2,moms3);
}

EvtComplex KspipiModel::amp_Belle2010(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3){
  EvtComplex amp(-2.537,0.923);
  SetMomenta(0,p4_p,moms1,moms2);//K*(892)
  SetMomenta(1,p4_p,moms1,moms2);//K0*(1430)
  SetMomenta(2,p4_p,moms1,moms2);//K2*(1430)
  SetMomenta(3,p4_p,moms1,moms2);//K*(1680)
  SetMomenta(4,p4_p,moms1,moms2);//K*(1410)

  SetMomenta(5,p4_p,moms1,moms3);//DCS K*(892)
  SetMomenta(6,p4_p,moms1,moms3);//DCS K0*(1430)
  SetMomenta(7,p4_p,moms1,moms3);//DCS K2*(1430)
  SetMomenta(8,p4_p,moms1,moms3);//DCS K*(1680)
  SetMomenta(9,p4_p,moms1,moms3);//DCS K*(1410)

  SetMomenta(10,p4_p,moms3,moms2);//Rho
  SetMomenta(11,p4_p,moms3,moms2);//Omega
  SetMomenta(12,p4_p,moms3,moms2);//f0(980)
  SetMomenta(13,p4_p,moms3,moms2);//f0(1370)
  SetMomenta(14,p4_p,moms3,moms2);//f2(1270)
  SetMomenta(15,p4_p,moms3,moms2);//Rho(1450)
  SetMomenta(16,p4_p,moms3,moms2);//sigma1
  SetMomenta(17,p4_p,moms3,moms2);//sigma2
  for(int i=0; i<ResNum(); i++){ amp += const_cast<EvtResonance2*>(Res(i))->resAmpl();}
  return amp;
}
