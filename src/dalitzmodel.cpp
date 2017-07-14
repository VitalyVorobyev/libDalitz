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

using compld = std::complex<double>;
using vectcd = std::vector<compld>;

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;

DalitzModel::DalitzModel(double mm, double ma, double mb, double mc) :
    AbsDalitzModel(mm, ma, mb, mc), use_subset(false) {
    UnitAmps(true);
}

bool DalitzModel::IsVetoed(double mAB, double mAC) const {
    for (auto i = 0u; i < m_veto_v.size(); i++) {
        if ((*m_veto_v[i])(mAB, mAC)) return true;
    }
    return false;
}

void DalitzModel::InitNames() {
    vector<string> res_names;
    for (auto* res : m_res_v)
        res_names.push_back(res->Name());
    SetResNames(res_names);
}

auto DalitzModel::GetResAmp(const DalitzPlotObject* res,
                            double mAB, double mAC) const {
    switch (res->Path()) {
    case DalitzResonance::AB: return res->evaluate(mAC, m3sq(mAC, mAB));
    case DalitzResonance::AC: return res->evaluate(mAB, m3sq(mAC, mAB));
    case DalitzResonance::BC: return res->evaluate(mAC, mAB);
    }
    cerr << "DalitzModel::GetResAmp: wrong path " << res->Path() << endl;
    return compld(0., 0.);
}

void DalitzModel::GetResVals(vectcd* resv, double mABsq, double mACsq) const {
    resv->clear();
    if (IsVetoed(mABsq, mACsq)) {
        resv->resize(m_res_v.size(), 0);
        return;
    }
    if (use_subset) {  // If some subset is chosen
        for (auto res : m_rlist)
            resv->push_back(GetResAmp(m_res_v[res], mABsq, mACsq));
    } else {
        for (auto res : m_res_v)
            resv->push_back(GetResAmp(res, mABsq, mACsq));
    }
}

auto DalitzModel::ResAmp(uint32_t n, double mAB, double mAC) const {
    if (n < 0 || n >= m_res_v.size()) return compld(0);
    return GetResAmp(m_res_v[n], mAB, mAC);
}

compld DalitzModel::GetResVal(double mABsq, double mACsq, int resnum) const {
    return ResAmp(resnum, mABsq, mACsq);
}

bool DalitzModel::SetRVec(const vector<uint32_t>& rlist) {
    m_rlist.clear(); use_subset = false;
    cout << "Resonances subset: ";
    for (auto i = 0u; i < rlist.size(); i++) {
        if (rlist[i] < m_res_v.size()) {
            cout << rlist[i] << " " << m_res_v[rlist[i]]->Name() << ", ";
            const uint32_t index = rlist[i];
            m_rlist.push_back(index);
        } else {
            cerr << endl << "Wrong resonance index " << rlist[i]
                 << " (should be less than " << m_res_v.size() << ")"<< endl;
        }
    }
    cout << endl;
    use_subset = true;
    return true;
}

int DalitzModel::Norm(double *val, double *err, uint64_t nc) {
    DalitzMCIntegral mcint(this);
    return mcint.CalcIntegral(val, err, nc);
}
