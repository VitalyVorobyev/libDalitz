/** Copyright 2017 Vitaly Vorobyev
 ** @file abssymdalitzmodel.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/abssymdalitzmodel.h"

#include <iostream>
#include <complex>
#include <fstream>
#include <iomanip>
#include <algorithm>

typedef AbsSymDalitzModel ASModel;
typedef std::complex<double> compld;
typedef std::string str;

using std::cerr;
using std::cout;
using std::endl;
using std::arg;
using std::max;

const int sign = -1;
const double ASModel::m_2pi = 2.*M_PI;

ASModel::AbsSymDalitzModel(const double& mm, const double& ma, const double& mb,
                 const double& delmin, const double& delmax) :
    AbsDalitzModel(mm, ma, mb, mb),
    del_min(delmin), del_max(delmax), m_nbins(8) {}

int ASModel::sign(const double& x) {
    return x >= 0 ? 1 : -1;
}

int ASModel::bin(const double& mAB, const double& mAC) const {
    if (!IsInPlot(mAB, mAC)) return 0;
    return bin(mAB, mAC, delta(mAB, mAC));
}

int ASModel::bin(const double& mAB, const double& mAC,
                 const double& dphi) const {
    return (mAB < mAC) ? bin(dphi, false) : -bin(dphi, true);
}

int ASModel::bin(const double& dphi, const bool sign) const {
    double del = sign ? - dphi : dphi;
    if (del < del_min) del += m_2pi;
    if (del > del_max) del -= m_2pi;
    return (static_cast<int>(m_nbins * del / m_2pi + 0.5) % m_nbins) + 1;
}

void ASModel::ppdelt(const double& mp, const double& mm,
                     double* const pp, double* const pn,
                     double* const del) const {
    const compld ap = Amp(mp, mm);
    const compld an = Amp(mm, mp);
    *del = delta(ap, an);
    *pp = norm(ap);
    *pn = norm(an);
}

double ASModel::delta(const double& mp, const double& mm) const {
    return delta(arg(Amp(mp, mm)), arg(Amp(mm, mp)));
}

double ASModel::delta(const compld& ap, const compld& an) const {
    double del = arg(ap) - arg(an);
    if (del < del_min) return del + m_2pi;
    if (del > del_max) return del - m_2pi;
    return del;
}

void ASModel::TabulateSymABAC(const str& fname,
                              const unsigned grid_size) const {
    std::ofstream file(fname, std::ofstream::out);
    if (!file.is_open()) {
        cerr << "Can't open file " << fname << endl;
        return;
    }
    file << "mAB x mAC, GridSize " << grid_size << " (ASModel)"
         << endl << AsText();
    file << std::setprecision(7) << std::fixed;
    // dmAB == dmAC because of symmetry
    const double dm = (mABsq_max() - mABsq_min()) / grid_size;
    for (double mAB = mABsq_min(); mAB < mABsq_max(); mAB += dm) {
        double mACmax, mACmin;
        mACsqRange_AB(mAB, &mACmin, &mACmax);
        for (double mAC = max(mAB, mACmin) ; mAC < mACmax; mAC += dm) {
            const compld ap = Amp(mAB, mAC);
            const compld an = Amp(mAC, mAB);
            if (std::isnan(std::real(ap)) || std::isnan(std::real(an))) {
                cerr << "TabulateSymABAC: nan detected: "
                     << mAB << " " << mAC << " "
                     << ap << " " << an << endl;
                continue;
            }
            file << mAB << " " << mAC << " "
                 << ap << " " << an << " " << bin(mAB, mAC, delta(ap, an)) << endl;
        }
    }
    file.close();
}

void ASModel::TabulateSymABBC(const str& fname,
                              const unsigned grid_size) const {
    std::ofstream file(fname, std::ofstream::out);
    if (!file.is_open()) {
        cerr << "Can't open file " << fname << endl;
        return;
    }
    file << "mAB x mBC, GridSize " << grid_size << " (SymModel)"
         << endl << AsText();
    file << std::setprecision(7) << std::fixed;
    const double dmAB = (mABsq_max() - mABsq_min()) / grid_size;
    const double dmBC = (mBCsq_max() - mBCsq_min()) / grid_size;
    for (double mAB = mABsq_min(); mAB < mABsq_max(); mAB += dmAB) {
        double mBCmax, mBCmin;
        mBCsqRange_AB(mAB, &mBCmin, &mBCmax);
        for (double mBC = mBCmin; mBC < mBCmax; mBC += dmBC) {
            const double mAC = m3sq(mAB, mBC);
            const compld ap = Amp(mAB, mAC);
            const compld an = Amp(mAC, mAB);
            if (std::isnan(std::real(ap)) || std::isnan(std::real(an))) {
                cerr << "TabulateSymABBC: nan detected: "
                     << mAB << " " << mAC << " " << mBC << " "
                     << ap << " " << an << endl;
                continue;
            }
            file << mAB << " " << mBC << " "
                 << ap << " " << an << " " << bin(mAB, mAC, delta(ap, an)) << endl;
        }
    }
    file.close();
}
