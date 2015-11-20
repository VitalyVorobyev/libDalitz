#include "dalitzresonance.h"

#include "dalitzphasespace.h"
#include "relbreitwigner.h"
#include "blattweisskopf.h"
#include "gounarissakurai.h"
#include "nrpropagator.h"
#include "buggpropagator.h"
#include "virtualdstarpropagator.h"

#include <math.h>

DalitzResonance::DalitzResonance(const std::string& name,const int PropType, const int WidthType,const double& mmo, const double& mca, const double& mcb, const double& mcc, const double& mres, const double& wres, const int spin, const EvtComplex& camp):
  DalitzPlotObject(name,camp),
  m_ang_amp(new ResDecayAngularDistribution(spin,mmo,mca,mcb,mcc,mres)),
  m_ptype(PropType),m_wtype(WidthType)
// Constuctor for RBW and GS resonances
{
  const ResDecayAngularDistribution& angamp = *m_ang_amp;
  const double pResRef = DalitzPhaseSpace::pResSq(angamp.mResSq(),angamp.mChASq(),angamp.mChBSq());
  switch(PropType){
  case ResPropType::RBW:
    if(WidthType == VarWType::Flatte){
      m_prop = new RelBreitWigner(wres,mres,pResRef,spin,VarWType::Flatte);
    } else{
      m_prop = new RelBreitWigner(wres,mres,pResRef,spin,WidthType == VarWType::Const ? true : false);
    }
    break;
  case ResPropType::GS:
    m_prop = new GounarisSakurai(wres,mres,pResRef,WidthType == VarWType::Const ? true : false);
    break;
  case ResPropType::Flatte:
    m_prop = new RelBreitWigner(wres,mres,pResRef,spin,VarWType::Flatte);
    break;
  default:
    std::cout << "DalitzResonance: wrong proparator type for " << name;
    std::cout << ". Should be " << ResPropType::RBW << " or " << ResPropType::GS << std::endl;
    return;
    break;
  }
  const double pMoRef = DalitzPhaseSpace::pResSq(angamp.mMotSq(),angamp.mResSq(),angamp.mChCSq());
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

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& alpha, const EvtComplex& camp):
  DalitzPlotObject(name,camp),
  m_ang_amp(new ResDecayAngularDistribution(0,mmo,mca,mcb,mcc,0)),
  m_ptype(PropType)
{
  if(PropType != ResPropType::NR && PropType != ResPropType::Flatte){
    std::cout << "DalitzResonance: wrong propagator type specified in NR of Flatte constructor for " << name << std::endl;
    return;
  }
  m_wtype = PropType == ResPropType::Flatte ? VarWType::Flatte : VarWType::Const;
  m_prop  = PropType == ResPropType::Flatte ? (AbsPropagator*) new RelBreitWigner(1,alpha,1,0,VarWType::Flatte) : new NRPropagator(alpha);
  const ResDecayAngularDistribution& angamp = *m_ang_amp;
  const double pResRef = DalitzPhaseSpace::pResSq(angamp.mResSq(),angamp.mChCSq(),angamp.mChBSq());
  const double pMoRef  = DalitzPhaseSpace::pResSq(angamp.mMotSq(),angamp.mResSq(),angamp.mChBSq());
  m_mff = new BlattWeisskopf(0,1.6,pMoRef);
  m_rff = new BlattWeisskopf(0,1.6,pResRef);
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& alpha,const double& amp,const double& phi):
  DalitzResonance(name,PropType,mmo,mca,mcb,mcc,alpha,amp*EvtComplex(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const EvtComplex& camp):
  DalitzPlotObject(name,camp),
  m_ang_amp(new ResDecayAngularDistribution(0,mmo,mca,mcb,mcc,0.5)),
  m_prop(new BuggPropagator()),
  m_ptype(PropType),m_wtype(VarWType::Const)
{
  if(PropType != ResPropType::Bugg){
    std::cout << "DalitzResonance: wrong propagator type specified in Bugg constructor for " << name << std::endl;
    return;
  }
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& amp,const double& phi)
  : DalitzResonance(name,PropType,mmo,mca,mcb,mcc,amp*EvtComplex(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& beta1,const double& beta2, const EvtComplex& camp):
  DalitzPlotObject(name,camp),
  m_ang_amp(new ResDecayAngularDistribution(1,mmo,mca,mcb,mcc,2.01)),
  m_prop(new VirtualDstarPropagator(beta1,beta2)),m_ptype(PropType)
{
  if(PropType != ResPropType::VDst){
    std::cout << "DalitzResonance: wrong propagator type specified in VDst constructor for " << name << std::endl;
    return;
  }
  const ResDecayAngularDistribution& angamp = *m_ang_amp;
  const double pResRef = DalitzPhaseSpace::pResSq(angamp.mResSq(),angamp.mChASq(),angamp.mChBSq());
  const double pMoRef  = DalitzPhaseSpace::pResSq(angamp.mMotSq(),angamp.mResSq(),angamp.mChCSq());
  m_mff = new BlattWeisskopf(0,1.6,pMoRef);
  m_rff = new BlattWeisskopf(0,1.6,pResRef);
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const double& mmo, const double& mca, const double& mcb, const double& mcc,const double& beta1,const double& beta2, const double& amp,const double& phi):
  DalitzResonance(name,PropType,mmo,mca,mcb,mcc,beta1,beta2,amp*EvtComplex(cos(phi),sin(phi)))
{
}

EvtComplex DalitzResonance::evaluate(const double& mACsq,const double& mBCsq){
  const ResDecayAngularDistribution& angamp = *m_ang_amp;
  const double mABsq  = DalitzPhaseSpace::mBC(angamp.mMotSq(),angamp.mChASq(),angamp.mChBSq(),angamp.mChCSq(),mACsq,mBCsq);
  const double pMoAB  = DalitzPhaseSpace::pResSq(angamp.mMotSq(),mABsq,angamp.mChCSq());
  const double pResAB = DalitzPhaseSpace::pResSq(mABsq,angamp.mChASq(),angamp.mChBSq());
  return CAmp()*(*m_mff)(pMoAB)*(*m_rff)(pResAB)*angamp(mACsq,mBCsq)*(*m_prop)(mABsq,pResAB);
}
