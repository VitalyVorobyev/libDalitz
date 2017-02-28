/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzmodel.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/dalitzmodel.h"

#include <iostream>

#include "./dalitzmcintegral.h"
#include "./dalitzresonance.h"

typedef std::complex<double> compld;
typedef std::vector<compld> vectcd;

using std::cout;
using std::endl;

DalitzModel::DalitzModel(const double& mm, const double& ma,
                         const double& mb, const double& mc) :
    AbsDalitzModel(mm, ma, mb, mc), use_subset(false) {}

bool DalitzModel::IsVetoed(const double& mAB, const double& mAC) const {
    for (unsigned i=0; i < m_veto_v.size(); i++) {
        if ((*m_veto_v[i])(mAB, mAC)) return true;
    }
    return false;
}

compld DalitzModel::Amp(const double& mAB, const double& mAC) const {
    compld amp(0., 0.);
    if (IsVetoed(mAB, mAC)) return amp;
    if (use_subset) {  // If some subset is chosen
        for (unsigned i=0; i < m_rlist.size(); i++)
            amp += GetResAmp(m_res_v[m_rlist[i]], mAB, mAC);
    } else {
        for (unsigned i=0; i < m_res_v.size(); i++)
            amp += GetResAmp(m_res_v[i], mAB, mAC);
    }
    return amp;
}

void DalitzModel::GetResVals(vectcd* resv, const double& mABsq,
                             const double& mACsq) const {
    resv->clear();
    if (IsVetoed(mABsq, mACsq)) {
        resv->resize(m_res_v.size(), 0);
        return;
    }
    if (use_subset) {  // If some subset is chosen
        for (unsigned i=0; i < m_rlist.size(); i++)
            resv->push_back(GetResAmp(m_res_v[m_rlist[i]], mABsq, mACsq));
    } else {
        for (unsigned i=0; i < m_res_v.size(); i++)
            resv->push_back(GetResAmp(m_res_v[i], mABsq, mACsq));
    }
}

compld DalitzModel::ResAmp(const unsigned n, const double& mAB,
                           const double& mAC) const {
    if (n < 0 || n >= m_res_v.size()) return 0;
    return GetResAmp(m_res_v[n], mAB, mAC);
}

compld DalitzModel::GetResAmp(const DalitzPlotObject* res, const double& mAB,
                              const double& mAC) const {
    switch (res->Path()) {
    case DalitzResonance::AB: return res->evaluate(mAC, m3sq(mAC, mAB));
    case DalitzResonance::AC: return res->evaluate(mAB, m3sq(mAC, mAB));
    case DalitzResonance::BC: return res->evaluate(mAC, mAB);
    }
    return compld(0., 0.);
}

compld DalitzModel::GetAmplitudes(vectcd* vec, const double& mAB,
                                  const double& mAC) const {
    vec->clear();
    compld amp(0., 0.);
    for (unsigned i=0; i < m_res_v.size(); i++) {
        vec->push_back(GetResAmp(m_res_v[i], mAB, mAC));
        amp += vec->at(i);
    }
    return amp;
}

bool DalitzModel::SetRVec(const std::vector<unsigned>& rlist) {
    m_rlist.clear(); use_subset = false;
    cout << "Resonances subset: ";
    for (unsigned i=0; i < rlist.size(); i++) {
        if (rlist[i] < m_res_v.size()) {
            cout << rlist[i] << " " << m_res_v[rlist[i]]->Name() << ", ";
            const unsigned index = rlist[i];
            m_rlist.push_back(index);
        } else {
            cout << endl << "Wrong resonance index " << rlist[i]
                 << " (should be less than " << m_res_v.size() << ")"<< endl;
        }
    }
    cout << endl;
    use_subset = true;
    return true;
}

int DalitzModel::Norm(double *val, double *err, const uint64_t nc) {
    DalitzMCIntegral mcint(this);
    return mcint.CalcIntegral(val, err, nc);
}

