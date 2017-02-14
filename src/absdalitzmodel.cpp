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

typedef AbsDalitzModel AbsDM;
typedef std::complex<double> compld;
typedef std::vector<double> vectd;
typedef std::vector<int> vecti;
typedef std::vector<compld> vectcd;
typedef std::string str;

using std::cout;
using std::endl;
using std::fabs;
using std::to_string;
using std::abs;
using std::arg;
using std::norm;
using std::real;

AbsDalitzModel::AbsDalitzModel(const DalitzPhaseSpace& phsp) :
    DalitzPhaseSpace(phsp) {}

AbsDalitzModel::AbsDalitzModel(const double& mM, const double& mA,
                               const double& mB, const double& mC) :
    DalitzPhaseSpace(mM, mA, mB, mC) {}

compld AbsDM::Amp(const double& mABsq, const double& mACsq) const {
    vectcd resvals;
    compld amp = 0;
    GetAmpVals(&resvals, mABsq, mACsq);
    for (unsigned i=0; i < m_ampl.size(); i++) amp += m_ampl[i]*resvals[i];
    return amp;
}

std::vector<MnPar> AbsDM::MnPars(void) const {
    const double phmin = -2.0*M_PI;
    const double phmax =  2.0*M_PI;
    const double pherr =  0.1*M_PI;
    std::vector<MnPar> parv;
    if (m_parstate.size() != 2*m_amp_names.size()) {
        cout << "Wrong size of m_parstate: " << m_parstate.size()
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

double AbsDM::P(const double& mABsq, const double& mACsq) const {
    return norm(Amp(mABsq, mACsq));
}
double AbsDM::Arg(const double& mABsq, const double& mACsq) const {
    return arg(Amp(mABsq, mACsq));
}

compld AbsDM::GetAmplitudes(vectcd* resv,
                            const double& mABsq, const double& mACsq) const {
    resv->clear();
    vectcd resvals;
    GetResVals(&resvals, mABsq, mACsq);
    if (m_ampl.size() == resvals.size()) {
        for (unsigned i=0; i < m_amp_signature.size(); i++)
            resv->push_back(m_ampl[i]*resvals[i]);
    } else {
        int res = 0;
        for (unsigned i=0; i < m_amp_signature.size(); i++) {
            for (unsigned j=0; j < m_amp_signature[i]; j++) {
                resv->push_back(m_ampl[i]*resvals[res++]);
            }
        }
    }
    compld sum(0);
    for (compld amp : *resv) sum += amp;
    return sum;
}

void AbsDM::GetAmpVals(vectcd* resv,
                       const double& mABsq, const double& mACsq) const {
    resv->clear();
    vectcd* resvals = new vectcd();
    GetResVals(resvals, mABsq, mACsq);
    if (m_ampl.size() == resvals->size()) {
        resv = resvals;
        return;
    } else {
        int resnum = 0;
        for (unsigned i=0; i < m_amp_signature.size(); i++) {
            resv->push_back(0);
            for (unsigned j=0; j < m_amp_signature[i]; j++) {
                resv->at(i) += resvals->at(resnum++);
            }
        }
    }
}

int AbsDM::OpenCachedIntegrals(const str& fname, const bool silent) {
    std::ifstream ifile(fname, std::ifstream::in);
    if (!ifile.is_open()) {
        cout << "Can't open file " << fname << endl;
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
    for (compld amp : m_ampl) out << abs(amp) << "*exp(i*" << arg(amp) << ") ";
    return out.str();
}

void AbsDM::ShowAmpls(void) const {
    cout << GetAmpStr() << endl;
}

double AbsDalitzModel::NormWithCache(void) const {
    if (!m_ampl.size() || m_ampl.size() != m_res_int.size()) {
        cout << "Can't calculate normalizations with cached integrals: "
             << m_ampl.size() << ", " << m_res_int.size() << endl;
        return 1;
    }
    double res = 0;
    for (unsigned i=0; i < m_ampl.size(); i++) {
        // I += |a_i|^2 * I_i
        res += norm(m_ampl[i])*real(m_res_int[i][i]);
        for (unsigned j=i+1; j < m_ampl.size(); j++) {
            // I += 2*Re( a_i * a_j^* * I_{ij} )
            res += 2.*real(m_ampl[i]*conj(m_ampl[j])*m_res_int[i][j]);
        }
    }
    return res;
}

void AbsDalitzModel::GetCoefficients(vectcd* coefv) const {*coefv = m_ampl;}

void AbsDalitzModel::SetResAreas(const vectd& ledge, const vectd& redge,
                                 const vecti& types) {
    m_res_areas.clear();
    for (unsigned i=0; i < types.size(); i++)
        m_res_areas.push_back(new DStrip(ledge[i], redge[i], types[i]));
}
