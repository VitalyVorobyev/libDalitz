#include "kspipimodel.h"
#include "consts.h"
#include "blattweisskopf.h"
#include "absvarwidth.h"

KspipiModel::KspipiModel(void):
KspipiModel(m_D0_Mass,m_Ks0_Mass,m_PI_Mass)
{
}
KspipiModel::KspipiModel(const double &md, const double &mks, const double& mpi):
SymDalitzModel(md,mks,mpi,-M_PI/8.,15.*M_PI/8.)
{
  SetABaxis("m_{+}^{2}, GeV^{2}/c^{4}");
  SetACaxis("m_{-}^{2}, GeV^{2}/c^{4}");
  SetBCaxis("m_{#pi#pi}^{2}, GeV^{2}/c^{4}");
  // ** A. Poluektov et al. Phys. Rev. D 81, 112002 – Published 16 June 2010 **
  const double dtr = 1./EvtConst::radToDegrees;
  BlattWeisskopf::m_r_meson     = 5.0;
  BlattWeisskopf::m_r_resonance = 1.5;

  // CF //
  std::cout << "VarWType::BW = " << VarWType::BW << std::endl;
  AddRes(new DalitzResonance("K*(892)",      ResPropType::RBW,VarWType::BW,this,ResPath::AB,0.8937,0.0484,1,1.638,133.2*dtr));
//  AddRes(new DalitzResonance("K0*(1430)",    ResPropType::RBW,this,ResPath::AB,1.4120,0.2940,0,2.210,358.9*dtr));
//  AddRes(new DalitzResonance("K2*(1430)",    ResPropType::RBW,this,ResPath::AB,1.4256,0.0985,2,0.890,314.8*dtr));
//  AddRes(new DalitzResonance("K*(1680)",     ResPropType::RBW,this,ResPath::AB,1.7170,0.3220,1,0.880, 82.0*dtr));
//  AddRes(new DalitzResonance("K*(1410)",     ResPropType::RBW,this,ResPath::AB,1.4140,0.2320,1,0.650,120.0*dtr));
//  // DCS //
//  AddRes(new DalitzResonance("K*(892) DCS",  ResPropType::RBW,this,ResPath::AC,.89166,0.0508,1,0.149,325.4*dtr));
//  AddRes(new DalitzResonance("K0*(1430) DCS",ResPropType::RBW,this,ResPath::AC,1.4120,0.2940,0,0.360, 87.0*dtr));
//  AddRes(new DalitzResonance("K2*(1430) DCS",ResPropType::RBW,this,ResPath::AC,1.4256,0.0985,2,0.230,275.0*dtr));
//  AddRes(new DalitzResonance("K*(1680) DCS", ResPropType::RBW,this,ResPath::AC,1.7170,0.3220,1,2.100,130.0*dtr));
//  AddRes(new DalitzResonance("K*(1410) DCS", ResPropType::RBW,this,ResPath::AC,1.4140,0.2320,1,0.420,253.0*dtr));
//  // CP
//  AddRes(new DalitzResonance("rho(770)",     ResPropType::RBW,this,ResPath::BC,0.7717,0.1490,1,1.000,0.000));
//  AddRes(new DalitzResonance("omega(782)",   ResPropType::RBW,this,ResPath::BC,.78265,.00849,1,.0343,112.0*dtr));
//  AddRes(new DalitzResonance("f0(980)",      ResPropType::RBW,this,ResPath::BC,0.9770,0.0500,0,0.385,207.3*dtr));
//  AddRes(new DalitzResonance("f0(1370)",     ResPropType::RBW,this,ResPath::BC,1.3100,0.2720,0,1.250, 69.0*dtr));
//  AddRes(new DalitzResonance("f2(1270)",     ResPropType::RBW,this,ResPath::BC,1.2754,0.1851,2,1.440,342.9*dtr));
//  AddRes(new DalitzResonance("rho(1450)",    ResPropType::RBW,this,ResPath::BC,1.4650,0.4000,1,0.490, 64.0*dtr));
//  AddRes(new DalitzResonance("sigma1",       ResPropType::RBW,this,ResPath::BC,0.5220,0.4530,0,1.560,214.0*dtr));
//  AddRes(new DalitzResonance("sigma2",       ResPropType::RBW,this,ResPath::BC,1.0330,0.0880,0,0.200,212.0*dtr));
//  // NR
//  AddRes(new DalitzResonance("NR",           ResPropType::NR,this,ResPath::BC,0,EvtComplex(-2.537,0.923)));
}

//EvtComplex KspipiModel::Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3){
//  return amp_Belle2010(p4_p,moms1,moms2,moms3);
//}

//EvtComplex KspipiModel::amp_Belle2010(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3){
//  EvtComplex amp(-2.537,0.923);
//  SetMomenta(0,p4_p,moms1,moms2);//K*(892)
//  SetMomenta(1,p4_p,moms1,moms2);//K0*(1430)
//  SetMomenta(2,p4_p,moms1,moms2);//K2*(1430)
//  SetMomenta(3,p4_p,moms1,moms2);//K*(1680)
//  SetMomenta(4,p4_p,moms1,moms2);//K*(1410)

//  SetMomenta(5,p4_p,moms1,moms3);//DCS K*(892)
//  SetMomenta(6,p4_p,moms1,moms3);//DCS K0*(1430)
//  SetMomenta(7,p4_p,moms1,moms3);//DCS K2*(1430)
//  SetMomenta(8,p4_p,moms1,moms3);//DCS K*(1680)
//  SetMomenta(9,p4_p,moms1,moms3);//DCS K*(1410)

//  SetMomenta(10,p4_p,moms3,moms2);//Rho
//  SetMomenta(11,p4_p,moms3,moms2);//Omega
//  SetMomenta(12,p4_p,moms3,moms2);//f0(980)
//  SetMomenta(13,p4_p,moms3,moms2);//f0(1370)
//  SetMomenta(14,p4_p,moms3,moms2);//f2(1270)
//  SetMomenta(15,p4_p,moms3,moms2);//Rho(1450)
//  SetMomenta(16,p4_p,moms3,moms2);//sigma1
//  SetMomenta(17,p4_p,moms3,moms2);//sigma2
//  for(int i=0; i<ResNum(); i++){ amp += const_cast<EvtResonance2*>(Res(i))->resAmpl();}
//  return amp;
//}
