#include "dalitzresonance.h"

#include "dalitzphasespace.h"
#include "relbreitwigner.h"
#include "blattweisskopf.h"
#include "virtualresff.h"
#include "gounarissakurai.h"
#include "nrpropagator.h"
#include "buggpropagator.h"
#include "virtualdstarpropagator.h"
//#include "virtualdstarpropagator2.h"
#include "rhoomegapropagator.h"

#include <cmath>
#include <iostream>

using namespace std;

typedef DalitzResonance DRes;

DalitzResonance::DalitzResonance(cstr& name,cint PropType, cint WidthType,cDPhSp* phsp, cint respath, cdouble& mres, cdouble& wres, cint spin, ccompld& camp):
  DalitzPlotObject(name,camp), m_ptype(PropType),m_wtype(WidthType),m_spin(spin),m_path(respath)
{
  cdouble pResRef = SetFFAngAmp(phsp,mres);
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
    cout << "DalitzResonance: wrong proparator type for " << name;
    cout << ". Should be " << ResPropType::RBW << " or " << ResPropType::GS << endl;
    return;
  }
  return;
}

DalitzResonance::DalitzResonance(cstr& name,cint PropType, cint WidthType,cDPhSp* phsp, cint respath, cdouble& mres, cdouble& wres, cint spin, cdouble& amp, cdouble& phi):
  DRes(name,PropType,WidthType,phsp,respath,mres,wres,spin,amp*compld(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(cstr& name,cint PropType,cDPhSp* phsp, cint respath, cdouble& mres, cdouble& wres, cint spin, ccompld& camp):
  DalitzResonance(name,PropType,VarWType::BW,phsp,respath,mres,wres,spin,camp)
{
}

DalitzResonance::DalitzResonance(cstr& name,cint PropType,cDPhSp* phsp, cint respath, cdouble& mres, cdouble& wres, cint spin, cdouble& amp, cdouble& phi):
  DalitzResonance(name,PropType,VarWType::BW,phsp,respath,mres,wres,spin,amp,phi)
{
}

DalitzResonance::DalitzResonance(cstr& name,cint PropType,cDPhSp* phsp, cint respath,cdouble& alpha, ccompld& camp):
  DalitzPlotObject(name,camp), m_ptype(PropType),m_spin(0),m_path(respath)
{
  if(PropType != ResPropType::NR && PropType != ResPropType::Flatte){
    cout << "DalitzResonance: wrong propagator type specified in NR of Flatte constructor for " << name << endl;
    return;
  }
  m_wtype = PropType == ResPropType::Flatte ? VarWType::Flatte : VarWType::Const;
  SetFFAngAmp(phsp,1.);
  m_prop  = PropType == ResPropType::Flatte ? (AbsPropagator*) new RelBreitWigner(1,alpha,1,0,VarWType::Flatte) : new NRPropagator(alpha);
}

DalitzResonance::DalitzResonance(cstr& name,cint PropType,cDPhSp* phsp, cint respath,cdouble& alpha,cdouble& amp,cdouble& phi):
  DalitzResonance(name,PropType,phsp,respath,alpha,amp*compld(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(cstr& name,cint PropType,cDPhSp* phsp, cint respath,const compld& camp):
  DalitzPlotObject(name,camp),m_ptype(PropType),m_path(respath)
{
  if(PropType != ResPropType::Bugg && PropType != ResPropType::VDst2){
    cout << "DalitzResonance: wrong propagator type specified in Bugg (of Dv*0(2010)) constructor for " << name << endl;
    return;
  }
  cdouble mres = PropType == ResPropType::VDst2 ? 2.01027 : 0.5;
  m_spin  = PropType == ResPropType::VDst2 ? 1 : 0;
  cdouble pResRef = SetFFAngAmp(phsp,mres,PropType == ResPropType::VDst2 ? true : false);
  m_wtype = PropType == ResPropType::VDst2 ? VarWType::BW : VarWType::Const;
  m_prop  = PropType == ResPropType::VDst2 ? (AbsPropagator*) new RelBreitWigner(0.0834,mres,pResRef,1,m_wtype) : (AbsPropagator*) new BuggPropagator();
}

DalitzResonance::DalitzResonance(cstr& name,cint PropType,cDPhSp* phsp, cint respath,cdouble& amp,cdouble& phi)
  : DalitzResonance(name,PropType,phsp,respath,amp*compld(cos(phi),sin(phi)))
{
}

DalitzResonance::DalitzResonance(cstr& name,cint PropType,cDPhSp* phsp, cint respath,cdouble& beta1,cdouble& beta2, ccompld& camp):
  DalitzPlotObject(name,camp),m_ptype(PropType),m_wtype(VarWType::Const),m_spin(1),m_path(respath)
{
  if(PropType != ResPropType::VDst && PropType != ResPropType::RhoOmega){
    cout << "DalitzResonance: wrong propagator type specified in VDst or RhoOmega constructor for " << name << endl;
    return;
  }
  SetFFAngAmp(phsp,PropType == ResPropType::VDst ? 2.01 : 0.77);
  m_prop = PropType == ResPropType::VDst ? (AbsPropagator*) new VirtualDstarPropagator(beta1,beta2) : (AbsPropagator*) new RhoOmegaPropagator(beta1,beta2);
}

DalitzResonance::DalitzResonance(cstr& name,cint PropType,cDPhSp* phsp, cint respath,cdouble& beta1,cdouble& beta2, cdouble& amp,cdouble& phi):
  DalitzResonance(name,PropType,phsp,respath,beta1,beta2,amp*compld(cos(phi),sin(phi)))
{
}

compld DalitzResonance::evaluate(cdouble& mACsq,cdouble& mBCsq) const{
  const ResDecayAngularDistribution& angamp = *m_ang_amp;
  cdouble mABsq = DalitzPhaseSpace::GetmBCsq(angamp.mMotSq(),angamp.mChASq(),angamp.mChBSq(),angamp.mChCSq(),mACsq,mBCsq);
  cdouble p     = DalitzPhaseSpace::pRes(angamp.mMotSq(),mABsq,angamp.mChCSq());
  cdouble q     = DalitzPhaseSpace::pRes(mABsq,angamp.mChASq(),angamp.mChBSq());
//  cdouble q = DalitzPhaseSpace::q(mABsq,m_ang_amp->mChA(),m_ang_amp->mChB());
//  cdouble p = DalitzPhaseSpace::p(mABsq,m_ang_amp->mMot(),m_ang_amp->mChC());
//  cout << Name() << " " << CAmp() << " " << (*m_mff)(p) << " " << (*m_rff)(q) << " " << angamp(mACsq,mBCsq) << " " << (*m_prop)(mABsq,q) << endl;
  return CAmp()*(*m_mff)(p)*(*m_rff)(q)*angamp(mACsq,mBCsq)*(*m_prop)(mABsq,q);
}

double DalitzResonance::SetFFAngAmp(cDPhSp* phsp, cdouble& mres, const bool vdst){
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
    cout << "DalitzResonance: wrong resonance path " << m_path << " for resonance " << Name() << endl;
    return 0;
  }

  cdouble q0 = DalitzPhaseSpace::pRes(m_ang_amp->mResSq(),m_ang_amp->mChASq(),m_ang_amp->mChBSq());
  cdouble p0 = DalitzPhaseSpace::pRes(m_ang_amp->mMotSq(),m_ang_amp->mResSq(),m_ang_amp->mChCSq());
//  cdouble q  = DalitzPhaseSpace::q(m_ang_amp->mResSq(),m_ang_amp->mChA(),m_ang_amp->mChB());
//  cdouble p  = DalitzPhaseSpace::p(m_ang_amp->mResSq(),m_ang_amp->mMot(),m_ang_amp->mChC());
//  cout << "q: " << q0 << " -> " << q << endl;
//  cout << "p: " << p0 << " -> " << p  << " " << m_ang_amp->mResSq() << " " << m_ang_amp->mMotSq() << " " << m_ang_amp->mChCSq() << endl;

//  cout << "SetFFAngAmp: " << sqrt(psqResRef) << " " << sqrt(psqMoRef) << " " << m_spin << endl;
//  cout << "  " << m_ang_amp->mMotSq() << " " << m_ang_amp->mResSq() << " " << m_ang_amp->mChASq() << " " << m_ang_amp->mChBSq() << " " << m_ang_amp->mChCSq() << endl;
  m_mff = new BlattWeisskopf(m_spin,p0,FFType::FFMeson);
  m_rff = !vdst ? (FormFactor*) new BlattWeisskopf(m_spin,q0,FFType::FFResonance) : (FormFactor*) new VirtualResFF(BlattWeisskopf::m_r_resonance,q0);
  return q0;
}
