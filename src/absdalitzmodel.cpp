/** Copyright 2017 Vitaly Vorobyev
 ** @file absdalitzmodel.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/absdalitzmodel.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <numeric>  // std::accumulate, std::inner_product

using AbsDM = AbsDalitzModel;
using compld = std::complex<double>;
using vectd = std::vector<double>;
using vecti = std::vector<int>;
using vectcd = std::vector<compld>;
using str = std::string;

using std::cout;
using std::cerr;
using std::endl;
using std::fabs;
using std::to_string;
using std::abs;
using std::arg;
using std::norm;
using std::real;

AbsDalitzModel::AbsDalitzModel(double mM, double mA,
                               double mB, double mC) :
    DalitzPhaseSpace(mM, mA, mB, mC), m_unit_amps(false) {}

compld AbsDM::Amp(double mABsq, double mACsq) const {
    vectcd resvals;
    GetAmpVals(&resvals, mABsq, mACsq);
    if (m_unit_amps) {
        return std::accumulate(resvals.begin(), resvals.end(), compld(0.));
    } else {
        return std::inner_product(m_ampl.begin(), m_ampl.end(), resvals.begin(), compld(0.));
    }
}

std::vector<MnPar> AbsDM::MnPars(void) const {
    const double phmin = -2.0*M_PI;
    const double phmax =  2.0*M_PI;
    const double pherr =  0.1*M_PI;
    std::vector<MnPar> parv;
    if (m_parstate.size() != 2*m_amp_names.size()) {
        cerr << "Wrong size of m_parstate: " << m_parstate.size()
             << ". Expected: " << 2*m_amp_names.size() << endl;
        return parv;
    }
    for (unsigned i=0; i < m_amp_names.size(); i++) {
        const double amp = fabs(m_ampl[i]);
        parv.push_back(MnPar("amp_" + m_amp_names[i], amp,
                             0.1*amp, amp*0.1, amp*10., m_parstate[2*i]));
        const double pha = arg(m_ampl[i]);
        parv.push_back(MnPar("arg_" + m_amp_names[i], pha,
                             pherr, phmin, phmax, m_parstate[2*i+1]));
    }
    return parv;
}

double AbsDM::P(double mABsq, double mACsq) const {
    return norm(Amp(mABsq, mACsq));
}

double AbsDM::Arg(double mABsq, double mACsq) const {
    return arg(Amp(mABsq, mACsq));
}

compld AbsDM::GetAmplitudes(vectcd* resv,
                            double mABsq, double mACsq) const {
    resv->clear();
    vectcd resvals;
    GetResVals(&resvals, mABsq, mACsq);
    if (m_ampl.size() == resvals.size()) {
        if (m_unit_amps) {
            *resv = resvals;
        } else {
            for (unsigned i=0; i < m_ampl.size(); i++)
                resv->push_back(m_ampl[i]*resvals[i]);
        }
    } else {
        int res = 0;
        for (unsigned i=0; i < m_amp_signature.size(); i++) {
            for (unsigned j=0; j < m_amp_signature[i]; j++) {
                resv->push_back(m_ampl[i]*resvals[res++]);
            }
        }
    }
    return std::accumulate(resv->begin(), resv->end(), compld(0.));
}

void AbsDM::GetAmpVals(vectcd* resv,
                       double mABsq, double mACsq) const {
    resv->clear();
    vectcd resvals;
    GetResVals(&resvals, mABsq, mACsq);
    if (m_ampl.size() == resvals.size()) {
        *resv = resvals;
        return;
    } else {
        int resnum = 0;
        for (unsigned i=0; i < m_amp_signature.size(); i++) {
            resv->push_back(0);
            for (unsigned j=0; j < m_amp_signature[i]; j++) {
                resv->at(i) += resvals[resnum++];
            }
        }
    }
}

int AbsDM::OpenCachedIntegrals(const str& fname, const bool silent) {
    std::ifstream ifile(fname, std::ifstream::in);
    if (!ifile.is_open()) {
        cerr << "Can't open file " << fname << endl;
        return -1;
    }
    if (!silent) cout << "Reading cashed normalization from file "
                      << fname << endl;
    str line;
    unsigned size = 0;
    getline(ifile, line);
    std::sscanf(line.c_str(), "Size: %u", &size);
    if (!silent) cout << size << " amplitudes" << endl;
    m_res_int.resize(size);
    for (unsigned i=0; i < size; i++) m_res_int[i].resize(size, 0);
    int num1, num2;
    double val_re, val_im, err_re, err_im;
    while (getline(ifile, line)) {
        std::sscanf(line.c_str(), "%d %d (%lf,%lf) +- (%lf,%lf)",
                    &num1, &num2, &val_re, &val_im, &err_re, &err_im);
        m_res_int[num1][num2] = compld(val_re, val_im);
    }
    if (!silent) {
        cout << "Success!" << endl;
        for (unsigned i=0; i < size; i++) {
            for (unsigned j=0; j < size; j++) {
                cout << i << " " << j << " " << m_res_int[i][j] << endl;
            }
        }
    }
    return 0;
}

str AbsDM::GetAmpStr(void) const {
    std::stringstream out;
    out.str(""); out << std::setprecision(2) << std::fixed;
    for (auto& amp : m_ampl) out << abs(amp) << "*exp(i*" << arg(amp) << ") ";
    return out.str();
}

str AbsDM::AsText(void) const {
    std::stringstream out; out.str("");
    out << "### " << m_title << " ###" << endl;
    out << mM() << " -> " << mA() << " " << mB() << " " << mC() << endl;
    out << "mAB: " << mABaxis << " (" << mABsq_min() << ", "
        << mABsq_max() << ")" << endl;
    out << "mAC: " << mACaxis << " (" << mACsq_min() << ", "
        << mACsq_max() << ")" << endl;
    out << "mBC: " << mBCaxis << " (" << mBCsq_min() << ", "
        << mBCsq_max() << ")" << endl;
    return out.str();
}

double AbsDM::NormWithCache(void) const {
    if (!m_ampl.size() || m_ampl.size() != m_res_int.size()) {
        cerr << "Can't calculate normalizations with cached integrals: "
             << m_ampl.size() << ", " << m_res_int.size() << endl;
        return 1;
    }
    double res = 0;
    for (auto i = 0u; i < m_ampl.size(); i++) {
        // I += |a_i|^2 * I_i
        res += norm(m_ampl[i])*real(m_res_int[i][i]);
        for (auto j = i + 1u; j < m_ampl.size(); j++) {
            // I += 2*Re( a_i * a_j^* * I_{ij} )
            res += 2.*real(m_ampl[i]*conj(m_ampl[j])*m_res_int[i][j]);
        }
    }
    return res;
}

const vectcd& AbsDM::GetCoefficients() const {return m_ampl;}

void AbsDM::SetResAreas(const vectd& ledge, const vectd& redge,
                        const vecti& types) {
    m_res_areas.clear();
    for (unsigned i=0; i < types.size(); i++)
        m_res_areas.emplace_back(new DStrip(ledge[i], redge[i], types[i]));
}

void AbsDM::TabulateABAC(const str& fname, const unsigned grid_size) const {
    std::ofstream file(fname, std::ofstream::out);
    if (!file.is_open()) {
        cerr << "Can't open file " << fname << endl;
        return;
    }
    file << "mAB x mAC, GridSize " << grid_size << endl << AsText();
    file << std::setprecision(7) << std::fixed;
    const double dmAB = (mABsq_max() - mABsq_min()) / grid_size;
    const double dmAC = (mACsq_max() - mACsq_min()) / grid_size;
    for (double mAB = mABsq_min(); mAB < mABsq_max(); mAB += dmAB) {
        double mACmax, mACmin;
        mACsqRange_AB(mAB, &mACmin, &mACmax);
        for (double mAC = mACmin; mAC < mACmax; mAC += dmAC) {
            const compld amp = Amp(mAB, mAC);
            if (std::isnan(std::real(amp))) {
                cerr << "TabulateABAC: nan detected: "
                     << mAB << " " << mAC << " " << mAC << amp << endl;
                continue;
            }
            file << mAB << " " << mAC << " " << amp << endl;
        }
    }
    file.close();
}

void AbsDM::TabulateABBC(const str& fname, const unsigned grid_size) const {
    std::ofstream file(fname, std::ofstream::out);
    if (!file.is_open()) {
        cerr << "Can't open file " << fname << endl;
        return;
    }
    file << "mAB x mBC, GridSize " << grid_size << endl << AsText();
    file << std::setprecision(7) << std::fixed;
    const double dmAB = (mABsq_max() - mABsq_min()) / grid_size;
    const double dmBC = (mBCsq_max() - mBCsq_min()) / grid_size;
    for (double mAB = mABsq_min(); mAB < mABsq_max(); mAB += dmAB) {
        double mBCmax, mBCmin;
        mBCsqRange_AB(mAB, &mBCmin, &mBCmax);
        for (double mBC = mBCmin; mBC < mBCmax; mBC += dmBC) {
            const compld amp = Amp(mAB, m3sq(mAB, mBC));
            if (std::isnan(std::real(amp))) {
                cerr << "TabulateABAC: nan detected: "
                     << mAB << " " << mBC << " " << amp << endl;
                continue;
            }
            file << mAB << " " << mBC << " "
                 << amp << endl;
        }
    }
    file.close();
}
