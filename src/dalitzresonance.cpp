/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzresonance.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/dalitzresonance.h"

#include <cmath>
#include <iostream>

#include "./dalitzphasespace.h"
#include "./relbreitwigner.h"
#include "./blattweisskopf.h"
#include "./virtualresff.h"
#include "./gounarissakurai.h"
#include "./nrpropagator.h"
#include "./buggpropagator.h"
#include "./virtualdstarpropagator.h"
#include "./rhoomegapropagator.h"

typedef DalitzResonance DRes;
typedef DalitzPhaseSpace DPhSp;
typedef std::string str;
typedef std::complex<double> compld;

using std::cout;
using std::endl;
using std::cos;
using std::sin;

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const int WidthType, const DPhSp* phsp,
                                 const int respath, const double& mres,
                                 const double& wres, const int spin,
                                 const compld& camp) :
    DalitzPlotObject(name, camp), m_ptype(PropType),
    m_wtype(WidthType), m_spin(spin), m_path(respath) {
    const double pResRef = SetFFAngAmp(phsp, mres);
    switch (PropType) {
    case ResPropType::RBW:
        m_prop = new RelBreitWigner(wres, mres, pResRef, spin, WidthType);
        break;
    case ResPropType::GS:
        m_prop = new GounarisSakurai(wres, mres, pResRef,
                     WidthType == VarWType::Const ? true : false);
        break;
    case ResPropType::Flatte:
        m_prop = new RelBreitWigner(wres, mres, pResRef, spin,
                                    VarWType::Flatte);
        break;
    default:
        cout << "DalitzResonance: wrong proparator type for " << name;
        cout << ". Should be " << ResPropType::RBW << " or "
             << ResPropType::GS << endl;
    }
}

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const int WidthType, const DPhSp* phsp,
                                 const int respath, const double& mres,
                                 const double& wres, const int spin,
                                 const double& amp, const double& phi) :
    DRes(name, PropType, WidthType, phsp, respath, mres, wres, spin,
         amp*compld(cos(phi), sin(phi))) {}

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const DPhSp* phsp, const int respath,
                                 const double& mres, const double& wres,
                                 const int spin, const compld& camp) :
    DalitzResonance(name, PropType, VarWType::BW, phsp, respath,
                    mres, wres, spin, camp) {}

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const DPhSp* phsp, const int respath,
                                 const double& mres, const double& wres,
                                 const int spin, const double& amp,
                                 const double& phi) :
    DalitzResonance(name, PropType, VarWType::BW, phsp, respath,
                    mres, wres, spin, amp, phi) {}

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const DPhSp* phsp, const int respath,
                                 const double& alpha, const compld& camp) :
    DalitzPlotObject(name, camp), m_ptype(PropType), m_spin(0),
    m_path(respath) {
    if (PropType != ResPropType::NR && PropType != ResPropType::Flatte) {
        cout << "DalitzResonance: wrong propagator type specified"
                "in NR of Flatte constructor for " << name << endl;
        return;
    }
    m_wtype = PropType == ResPropType::Flatte ? VarWType::Flatte :
                                                VarWType::Const;
    SetFFAngAmp(phsp, 1.);
    m_prop = PropType == ResPropType::Flatte ?
        reinterpret_cast<AbsPropagator*>(
             new RelBreitWigner(1, alpha, 1, 0, VarWType::Flatte)) :
        reinterpret_cast<AbsPropagator*>(new NRPropagator(alpha));
}

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const DPhSp* phsp, const int respath,
                                 const double& alpha, const double& amp,
                                 const double& phi) :
    DalitzResonance(name, PropType, phsp, respath, alpha,
                    amp*compld(cos(phi), sin(phi))) {}

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const DPhSp* phsp, const int respath,
                                 const compld& camp) :
    DalitzPlotObject(name, camp), m_ptype(PropType), m_path(respath) {
    if (PropType != ResPropType::Bugg && PropType != ResPropType::VDst2) {
        cout << "DalitzResonance: wrong propagator type specified"
                "in Bugg (of Dv*0(2010)) constructor for " << name << endl;
        return;
    }
    const double mres = PropType == ResPropType::VDst2 ? 2.01027 : 0.5;
    m_spin  = PropType == ResPropType::VDst2 ? 1 : 0;
    const double pResRef = SetFFAngAmp(phsp, mres,
                                PropType == ResPropType::VDst2 ? true : false);
    m_wtype = PropType == ResPropType::VDst2 ? VarWType::BW : VarWType::Const;
    m_prop  = PropType == ResPropType::VDst2 ?
        reinterpret_cast<AbsPropagator*>(
            new RelBreitWigner(0.0834, mres, pResRef, 1, m_wtype)) :
        reinterpret_cast<AbsPropagator*>(new BuggPropagator());
}

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const DPhSp* phsp, const int respath,
                                 const double& amp, const double& phi) :
    DalitzResonance(name, PropType, phsp, respath,
                    amp*compld(cos(phi), sin(phi))) {}

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const DPhSp* phsp, const int respath,
                                 const double& beta1, const double& beta2,
                                 const compld& camp) :
    DalitzPlotObject(name, camp), m_ptype(PropType),
    m_wtype(VarWType::Const), m_spin(1), m_path(respath) {
    if (PropType != ResPropType::VDst && PropType != ResPropType::RhoOmega) {
        cout << "DalitzResonance: wrong propagator type specified in"
                "VDst or RhoOmega constructor for " << name << endl;
        return;
    }
    SetFFAngAmp(phsp, PropType == ResPropType::VDst ? 2.01 : 0.77);
    m_prop = PropType == ResPropType::VDst ?
        reinterpret_cast<AbsPropagator*>(
            new VirtualDstarPropagator(beta1, beta2)) :
        reinterpret_cast<AbsPropagator*>(new RhoOmegaPropagator(beta1, beta2));
}

DalitzResonance::DalitzResonance(const str& name, const int PropType,
                                 const DPhSp* phsp, const int respath,
                                 const double& beta1, const double& beta2,
                                 const double& amp, const double& phi) :
    DalitzResonance(name, PropType, phsp, respath, beta1, beta2,
                    amp*compld(cos(phi), sin(phi))) {}

compld DalitzResonance::evaluate(const double& mACsq,
                                 const double& mBCsq) const {
    const ResDecayAngularDistribution& angamp = *m_ang_amp;
    const double mABsq =
            DalitzPhaseSpace::m3sq(angamp.mMotSq(), angamp.mChASq(),
                                   angamp.mChBSq(), angamp.mChCSq(),
                                   mACsq, mBCsq);
    const double p = DalitzPhaseSpace::pRes(angamp.mMotSq(),
                                            mABsq, angamp.mChCSq());
    const double q = DalitzPhaseSpace::pRes(mABsq, angamp.mChASq(),
                                            angamp.mChBSq());
    return CAmp()*(*m_mff)(p)*(*m_rff)(q)*angamp(mACsq, mBCsq)
            *(*m_prop)(mABsq, q);
}

double DalitzResonance::SetFFAngAmp(const DPhSp* phsp, const double& mres,
                                    const bool vdst) {
    switch (m_path) {
    case AB:
        m_ang_amp = new ResDecayAngularDistribution(m_spin,
                        phsp->mM(), phsp->mA(), phsp->mB(), phsp->mC(), mres);
        break;
    case AC:
        m_ang_amp = new ResDecayAngularDistribution(m_spin,
                        phsp->mM(), phsp->mA(), phsp->mC(), phsp->mB(), mres);
        break;
    case BC:
        m_ang_amp = new ResDecayAngularDistribution(m_spin,
                        phsp->mM(), phsp->mB(), phsp->mC(), phsp->mA(), mres);
        break;
    default:
        cout << "DalitzResonance: wrong resonance path " << m_path
             << " for resonance " << Name() << endl;
        return 0;
    }

    const double q0 = DalitzPhaseSpace::pRes(m_ang_amp->mResSq(),
                                             m_ang_amp->mChASq(),
                                             m_ang_amp->mChBSq());
    const double p0 = DalitzPhaseSpace::pRes(m_ang_amp->mMotSq(),
                                             m_ang_amp->mResSq(),
                                             m_ang_amp->mChCSq());
    m_mff = new BlattWeisskopf(m_spin, p0, BlattWeisskopf::FFMeson);
    m_rff = !vdst ?
        reinterpret_cast<FormFactor*>(
            new BlattWeisskopf(m_spin, q0, BlattWeisskopf::FFResonance)) :
        reinterpret_cast<FormFactor*>(
            new VirtualResFF(BlattWeisskopf::m_r_resonance, q0));
    return q0;
}
