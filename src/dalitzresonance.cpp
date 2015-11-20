#include "dalitzresonance.h"

#include "dalitzphasespace.h"
#include "relbreitwigner.h"
#include "blattweisskopf.h"
#include "gounarissakurai.h"
#include "nrpropagator.h"
#include "buggpropagator.h"
#include "virtualdstarpropagator.h"
#include "rhoomegapropagator.h"

#include <math.h>

DalitzResonance::DalitzResonance(const std::string& name,const int PropType, const int WidthType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const EvtComplex& camp):
  DalitzPlotObject(name,camp), m_ptype(PropType),m_wtype(WidthType),m_spin(spin),m_path(respath)
{
  const double pResRef = SetFFAngAmp(phsp,mres);
  switch(PropType){
  case ResPropType::RBW:
//    if(WidthType == VarWType::Flatte){
//      m_prop = new RelBreitWigner(wres,mres,pResRef,spin,VarWType::Flatte);
//    } else{
//      m_prop = new RelBreitWigner(wres,mres,pResRef,spin,WidthType == VarWType::Const ? true : false);
//    }
    m_prop = new RelBreitWigner(wres,mres,pResRef,spin,WidthType);
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
  }
  return;
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType, const int WidthType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const double& amp, const double& phi):
  DalitzResonance(name,PropType,WidthType,phsp,respath,mres,wres,spin,amp*EvtComplex(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const EvtComplex& camp):
  DalitzResonance(name,PropType,VarWType::BW,phsp,respath,mres,wres,spin,camp)
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath, const double& mres, const double& wres, const int spin, const double& amp, const double& phi):
  DalitzResonance(name,PropType,VarWType::BW,phsp,respath,mres,wres,spin,amp,phi)
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& alpha, const EvtComplex& camp):
  DalitzPlotObject(name,camp), m_ptype(PropType),m_spin(0),m_path(respath)
{
  if(PropType != ResPropType::NR && PropType != ResPropType::Flatte){
    std::cout << "DalitzResonance: wrong propagator type specified in NR of Flatte constructor for " << name << std::endl;
    return;
  }
  m_wtype = PropType == ResPropType::Flatte ? VarWType::Flatte : VarWType::Const;
  SetFFAngAmp(phsp,1.);
  m_prop  = PropType == ResPropType::Flatte ? (AbsPropagator*) new RelBreitWigner(1,alpha,1,0,VarWType::Flatte) : new NRPropagator(alpha);
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& alpha,const double& amp,const double& phi):
  DalitzResonance(name,PropType,phsp,respath,alpha,amp*EvtComplex(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const EvtComplex& camp):
  DalitzPlotObject(name,camp),m_prop(new BuggPropagator()),m_ptype(PropType),m_wtype(VarWType::Const),m_spin(0),m_path(respath)
{
  if(PropType != ResPropType::Bugg){
    std::cout << "DalitzResonance: wrong propagator type specified in Bugg constructor for " << name << std::endl;
    return;
  }
  SetFFAngAmp(phsp,0.5);
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& amp,const double& phi)
  : DalitzResonance(name,PropType,phsp,respath,amp*EvtComplex(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& beta1,const double& beta2, const EvtComplex& camp):
  DalitzPlotObject(name,camp),m_ptype(PropType),m_wtype(VarWType::Const),m_spin(1),m_path(respath)
{
  if(PropType != ResPropType::VDst && PropType != ResPropType::RhoOmega){
    std::cout << "DalitzResonance: wrong propagator type specified in VDst or RhoOmega constructor for " << name << std::endl;
    return;
  }
  SetFFAngAmp(phsp,PropType == ResPropType::VDst ? 2.01 : 0.77);
  m_prop = PropType == ResPropType::VDst ? (AbsPropagator*) new VirtualDstarPropagator(beta1,beta2) : (AbsPropagator*) new RhoOmegaPropagator(beta1,beta2);
}

DalitzResonance::DalitzResonance(const std::string& name,const int PropType,const DalitzPhaseSpace* phsp, const int respath,const double& beta1,const double& beta2, const double& amp,const double& phi):
  DalitzResonance(name,PropType,phsp,respath,beta1,beta2,amp*EvtComplex(cos(phi),sin(phi)))
{
}

EvtComplex DalitzResonance::evaluate(const double& mACsq,const double& mBCsq){
  const ResDecayAngularDistribution& angamp = *m_ang_amp;
  const double mABsq    = DalitzPhaseSpace::mBCsq(angamp.mMotSq(),angamp.mChASq(),angamp.mChBSq(),angamp.mChCSq(),mACsq,mBCsq);
  const double psqMoAB  = DalitzPhaseSpace::pResSq(angamp.mMotSq(),mABsq,angamp.mChCSq());
  const double psqResAB = DalitzPhaseSpace::pResSq(mABsq,angamp.mChASq(),angamp.mChBSq());
  return CAmp()*(*m_mff)(psqMoAB)*(*m_rff)(psqResAB)*angamp(mACsq,mBCsq)*(*m_prop)(mABsq,sqrt(psqResAB));
}

double DalitzResonance::SetFFAngAmp(const DalitzPhaseSpace *phsp, const double& mres){
  switch(m_path){
  case ResPath::AB:
    m_ang_amp = new ResDecayAngularDistribution(m_spin,phsp->mM(),phsp->mA(),phsp->mB(),phsp->mC(),mres);
    break;
  case ResPath::AC:
    m_ang_amp = new ResDecayAngularDistribution(m_spin,phsp->mM(),phsp->mA(),phsp->mC(),phsp->mB(),mres);
    break;
  case ResPath::BC:
    m_ang_amp = new ResDecayAngularDistribution(m_spin,phsp->mM(),phsp->mB(),phsp->mC(),phsp->mA(),mres);
    break;
  default:
    std::cout << "DalitzResonance: wrong resonance path " << m_path << " for resonance " << Name() << std::endl;
    return 0;
  }

  const double psqResRef = fabs(DalitzPhaseSpace::pResSq(m_ang_amp->mResSq(),m_ang_amp->mChASq(),m_ang_amp->mChBSq()));
  const double psqMoRef  = fabs(DalitzPhaseSpace::pResSq(m_ang_amp->mMotSq(),m_ang_amp->mResSq(),m_ang_amp->mChCSq()));
//  std::cout << "SetFFAngAmp: " << sqrt(psqResRef) << " " << sqrt(psqMoRef) << " " << m_spin << std::endl;
//  std::cout << "  " << m_ang_amp->mMotSq() << " " << m_ang_amp->mResSq() << " " << m_ang_amp->mChASq() << " " << m_ang_amp->mChBSq() << " " << m_ang_amp->mChCSq() << std::endl;
  m_mff = new BlattWeisskopf(m_spin,psqMoRef,FFType::FFMeson);
  m_rff = new BlattWeisskopf(m_spin,psqResRef,FFType::FFResonance);
  return sqrt(psqResRef);
}
