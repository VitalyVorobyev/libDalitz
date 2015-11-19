#include "dalitzresonance.h"

#include "dalitzphasespace.h"
#include "relbreitwigner.h"
#include "blattweisskopf.h"
#include "gounarissakurai.h"
#include "nrpropagator.h"
#include "buggpropagator.h"

#include <math.h>

const int ResPropType::NR       = 0;
const int ResPropType::RBW      = 1;
const int ResPropType::GS       = 2;
const int ResPropType::RhoOmega = 3;
const int ResPropType::Bugg     = 4;
const int ResPropType::VDst     = 5;

const int VarWType::Const  = 0;
const int VarWType::BW     = 1;
const int VarWType::GS     = 2;
const int VarWType::Flatte = 3;
const int VarWType::Bugg   = 4;

DalitzResonance::DalitzResonance(const std::string& name,const int PropType, const int WidthType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const EvtComplex& camp):
  DalitzPlotObject(name,camp),
  m_ang_amp(new ResDecayAngularDistribution(spin,mmo,mca,mcb,mcc,mres)),
  m_ptype(PropType),m_wtype(WidthType)
// Constuctor for RBW and GS resonances
{
  const ResDecayAngularDistribution& angamp = *m_ang_amp;
  const double pResRef = DalitzPhaseSpace::pResSq(angamp.mResSq(),angamp.mChCSq(),angamp.mChBSq());
  switch(PropType){
  case ResPropType::RBW:
    m_prop = new RelBreitWigner(wres,mres,pResRef,spin,WidthType == VarWType::Const ? true : false);
    break;
  case ResPropType::GS:
    m_prop = new GounarisSakurai(wres,mres,pResRef,WidthType == VarWType::Const ? true : false);
    break;
  default:
    std::cout << "DalitzResonance: wrong proparator type for " << name;
    std::cout << ". Should be " << ResPropType::RBW << " or " << ResPropType::GS << std::endl;
    return;
    break;
  }
  const double pMoRef = DalitzPhaseSpace::pResSq(angamp.mMotSq(),angamp.mResSq(),angamp.mChBSq());
  m_mff = new BlattWeisskopf(0,1.6,pMoRef);
  m_rff = new BlattWeisskopf(spin,1.6,pResRef);
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType, const int WidthType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const double& amp, const double& phi):
  DalitzResonance(name,PropType,WidthType,mmo,mca,mcb,mcc,mres,wres,spin,amp*EvtComplex(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const EvtComplex& camp):
  DalitzResonance(name,PropType,VarWType::BW,mmo,mca,mcb,mcc,mres,wres,spin,camp)
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const double& amp, const double& phi):
  DalitzResonance(name,PropType,VarWType::BW,mmo,mca,mcb,mcc,mres,wres,spin,amp,phi)
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& alpha, const EvtComplex& camp):
  DalitzPlotObject(name,camp),
  m_ang_amp(new ResDecayAngularDistribution(0,0,0,0,0,0)),
  m_ptype(PropType),m_wtype(VarWType::Const)
{
  if(PropType != ResPropType::NR){
    std::cout << "DalitzResonance: wrong propagator type specified in NR constructor for " << name << std::endl;
    return;
  }
  m_prop = new NRPropagator(alpha);
  const ResDecayAngularDistribution& angamp = *m_ang_amp;
  const double pResRef = DalitzPhaseSpace::pResSq(angamp.mResSq(),angamp.mChCSq(),angamp.mChBSq());
  const double pMoRef  = DalitzPhaseSpace::pResSq(angamp.mMotSq(),angamp.mResSq(),angamp.mChBSq());
  m_mff = new BlattWeisskopf(0,1.6,pMoRef);
  m_rff = new BlattWeisskopf(0,1.6,pResRef);
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& alpha,const double& amp,const double& phi):
  DalitzResonance(name,PropType,alpha,amp*EvtComplex(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const EvtComplex& camp):
  DalitzPlotObject(name,camp),
  m_ang_amp(new ResDecayAngularDistribution(0,0,0,0,0,0)),
  m_prop(new BuggPropagator()),
  m_ptype(PropType),m_wtype(VarWType::Const)
{
  if(PropType != ResPropType::Bugg){
    std::cout << "DalitzResonance: wrong propagator type specified in Bugg constructor for " << name << std::endl;
    return;
  }
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& amp,const double& phi)
  : DalitzResonance(name,PropType,amp*EvtComplex(cos(phi),sin(phi)))
{
}

//DalitzResonance::DalitzResonance(const std::string& name,const int PropType, const double& amp, const double& phi); /// Constuctor for NR and Bugg
//DalitzResonance::DalitzResonance(const std::string& name,const int PropType, const int WidthType, const double& mres, const EvtComplex& camp);               /// Constuctor for Flatte
//DalitzResonance::DalitzResonance(const std::string& name,const int PropType, const int WidthType, const double& mres, const double& amp, const double& phi); /// Constuctor for Flatte

EvtComplex DalitzResonance::evaluate(const double& mACsq,const double& mBCsq){
  const ResDecayAngularDistribution& angamp = *m_ang_amp;
  const double mABsq  = DalitzPhaseSpace::mBC(angamp.mMotSq(),angamp.mChASq(),angamp.mChBSq(),angamp.mChCSq(),mACsq,mBCsq);
  const double pMoAB  = DalitzPhaseSpace::pResSq(angamp.mMotSq(),mABsq,angamp.mChCSq());
  const double pResAB = DalitzPhaseSpace::pResSq(mABsq,angamp.mChCSq(),angamp.mChBSq());
  return CAmp()*(*m_mff)(pMoAB)*(*m_rff)(pResAB)*angamp(mACsq,mBCsq)*(*m_prop)(mABsq,pResAB);
}
