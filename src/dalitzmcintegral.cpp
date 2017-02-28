/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzmcintegral.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/dalitzmcintegral.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>  // max

typedef std::complex<double> compld;
typedef std::string str;
typedef std::vector<double> vectd;
typedef std::vector<compld> vectcd;

using std::cout;
using std::endl;
using std::sqrt;
using std::norm;
using std::pow;
using std::vector;
using std::ofstream;
using std::to_string;

const double DalitzMCIntegral::inv_sqrt_2pi = 0.3989422804014327;

DalitzMCIntegral::DalitzMCIntegral(const AbsDalitzModel *model) :
    RandomDalitzPoint(model), m_model(model), m_ncounts(1e6), m_prefix("") {}

double DalitzMCIntegral::GetDalitzPlotArea(const uint64_t& nc) const {
    const uint64_t Counts = nc ? nc : m_ncounts;
    const double mABsqMin = m_model->mABsq_min();
    const double mABsqMax = m_model->mABsq_max();

    const double mACsqMin = m_model->mACsq_min();
    const double mACsqMax = m_model->mACsq_max();

    uint64_t nhits = 0;
    double mABsq, mACsq;
    for (uint64_t count=0; count < Counts; count++) {
        GetUnconstrainedPoint(&mABsq, &mACsq);
        if (IsInPlot(mABsq, mACsq)) nhits++;
    }
    const double square = (mABsqMax - mABsqMin)*(mACsqMax - mACsqMin);
    const double p = static_cast<double>(nhits)/Counts;
    const double area  = square*p;
    const double darea = square*sqrt(p*(1.-p)/Counts);
    cout << "Dalitz Plot area = " << area << " +- " << darea
         << " (" << Counts << ") counts" << endl;
    return area;
}

int DalitzMCIntegral::CalcIntegral(double *val, double *err, const uint64_t& nc,
                                   const bool norm_flag) const {
    const uint64_t Counts = nc ? nc : m_ncounts;
    const double DPArea = m_model->Area() == 0 ? GetDalitzPlotArea(1e7) :
                                                 m_model->Area();
    *val = 0; *err = 0;
    double mABsq, mACsq;
    double res;
    cout << "DalitzMCIntegral: integrating with " << Counts
         << " points..." << endl;
    for (uint64_t i=0; i < Counts; i++) {
        if (!(i%100000)) cout << i << " counts val = " << *val
                              << " val2 = " << *err << endl;
        GetPoint(&mABsq, &mACsq);
        res = m_model->P(mABsq, mACsq);
        if (norm_flag) res /= m_model->Norm();
        *val += res/Counts;
        *err += res*res/Counts;
    }
    *err = sqrt((*err - pow(*val, 2))/Counts);
    *val *= DPArea;
    *err *= DPArea;
    cout << "Done! Int = " << *val << " +- " << *err
         << " (" << 100.*(*err)/(*val) << "%)" << endl;
    return 0;
}

void DalitzMCIntegral::CheckNormalization(const uint64_t &nc) const {
    cout << "Checking normalization... " << m_model->Norm() << endl;
    double val, err;
    CalcIntegral(&val, &err, nc, true);
}

int DalitzMCIntegral::CalcNormMap(str* ofile, const uint64_t& nc) const {
    const uint64_t Counts = nc ? nc : m_ncounts;
    const unsigned ncoef = m_model->AmpNum();
    const double DPArea  = m_model->Area() == 0 ? GetDalitzPlotArea(1e7) :
                                                  m_model->Area();
    // DPArea = 174.2 for B0 -> D0 pi+ pi-
    vector<vectcd> int_matrix;
    vector<vectcd> err_matrix;
    vector<vectd> err_matrix_real;
    vector<vectd> err_matrix_imag;
    int_matrix.resize(ncoef);
    err_matrix.resize(ncoef);
    err_matrix_real.resize(ncoef);
    err_matrix_imag.resize(ncoef);
    for (unsigned i=0; i < int_matrix.size(); i++) {
        err_matrix[i].resize(ncoef, 0);
        int_matrix[i].resize(ncoef, 0);
        err_matrix_real[i].resize(ncoef, 0);
        err_matrix_imag[i].resize(ncoef, 0);
    }
    vectcd res_amp_vec;
    double mABsq;
    double mACsq;
    for (uint64_t count=0; count < Counts; count++) {
        if (count && !(count%100000)) cout << count << " counts" << endl;
        GetPoint(&mABsq, &mACsq);
        m_model->GetAmpVals(&res_amp_vec, mABsq, mACsq);
        for (unsigned i=0; i < ncoef; i++) {
            for (unsigned j=i; j < ncoef; j++) {
                const compld cvalue = res_amp_vec[i]*conj(res_amp_vec[j]);
                int_matrix[i][j] += cvalue/compld(Counts);
                err_matrix_real[i][j] += pow(real(cvalue), 2)/Counts;
                err_matrix_imag[i][j] += pow(imag(cvalue), 2)/Counts;
            }
        }
    }
    for (unsigned i=0; i < ncoef; i++) {
        for (unsigned j=ncoef-1; j >= i; j--) {
            const double reerr = err_matrix_real[i][j];
            const double imerr = err_matrix_imag[i][j];
            const compld val      = int_matrix[i][j];
            int_matrix[i][j]     *= DPArea;
            err_matrix_real[i][j] =
                    DPArea*sqrt((reerr - pow(real(val), 2))/Counts);
            err_matrix_imag[i][j] =
                    DPArea*sqrt((imerr - pow(imag(val), 2))/Counts);
            err_matrix[i][j] =
                    compld(err_matrix_real[i][j], err_matrix_imag[i][j]);
            if (i != j) {
                int_matrix[j][i]      = int_matrix[i][j];
                err_matrix[j][i]      = err_matrix[i][j];
                err_matrix_real[j][i] = err_matrix_real[i][j];
                err_matrix_imag[j][i] = err_matrix_imag[i][j];
            }
        }
    }
    ShowIntMatrix(int_matrix, err_matrix);
    *ofile = "params/NaiveNormMap" + to_string(Counts);
    if (m_prefix != "") *ofile += "_" + m_prefix;
    *ofile += ".txt";
    WriteIntMatrix(*ofile, int_matrix, err_matrix);
    return 0;
}

double DalitzMCIntegral::CalcBranchings(vectd* brvec, vectd* brerr,
                                        const uint64_t& nc) const {
    const uint64_t Counts = nc ? nc : m_ncounts;
    const unsigned NRes = m_model->ResNum();
    brvec->resize(NRes, 0); brerr->resize(NRes, 0);
    vectcd amps;

    double mAB, mAC;
    double totint = 0;
    double toterr = 0;

    const double scale = 1.e5;
    cout << "DalitzMCIntegral: calculating branchings with "
         << Counts << " points..." << endl;
    for (uint64_t i=0; i < Counts; i++) {
        if (!(i%100000)) cout << i << " counts" << endl;
        GetPoint(&mAB, &mAC);
        const double res = scale*norm(m_model->GetAmplitudes(&amps, mAB, mAC));
        totint += res/Counts;
        toterr += res*res/Counts;
        for (unsigned j=0; j < NRes; j++) {
            const double res2 = scale*norm(amps[j]);
            brvec->at(j) += res2/Counts;
            brerr->at(j) += res2*res2/Counts;
        }
    }
    toterr = sqrt((toterr - totint*totint)/Counts);
    cout << "Full integral = " << totint << " +- " << toterr
         << " (" << 100.*toterr/totint << "%)" << endl;
    cout << std::setprecision(2) << std::scientific;
    for (unsigned j=0; j < NRes; j++) {
        brerr->at(j) = sqrt((brerr->at(j)-brvec->at(j)*brvec->at(j))/Counts);
        cout << " int = " << brvec->at(j) << " +- " << brerr->at(j);
        cout << " (" << (brvec->at(j) > 0 ? 100.*brerr->at(j)/brvec->at(j) : 0)
             << "%) for " << m_model->ResName(j) << endl;
    }
    double BrSum = 0;
    cout << endl;
    for (unsigned j=0; j < NRes; j++) {
        brvec->at(j) /= totint;
        const double err1 = brerr->at(j)/totint;
        const double err2 = brvec->at(j)*toterr/totint;
        brvec->at(j) *= 100;
        brerr->at(j) = sqrt(err1*err1+err2*err2)*100.;
        BrSum += brvec->at(j);
        cout << "Br = (" << brvec->at(j) << " +- " << brerr->at(j);
        cout << ")%, rel err " << (brvec->at(j) > 0 ?
             100.*brerr->at(j)/brvec->at(j) : 0);
        cout << "% for " << m_model->ResName(j) << endl;
    }
    cout << "Branchings sum: " << BrSum << endl;
    return totint;
}

vector<vectcd> DalitzMCIntegral::SmartNormMap(const uint64_t& nc) const {
    const uint64_t Counts = nc ? nc : m_ncounts;
    const unsigned ncoef  = m_model->AmpNum();

    vector<vectcd> int_matrix;
    int_matrix.resize(ncoef);
    for (unsigned i = 0; i < int_matrix.size(); i++)
        int_matrix[i].resize(ncoef, 0);

    for (unsigned i=0; i < ncoef; i++) {
        cout << "res1: " << i+1 << endl;
        int_matrix[i][i] = GetGaussIntegral(i, Counts);
        for (unsigned j=i+1; j < ncoef; j++) {
            cout << "res2: " << j+1 << endl;
            int_matrix[i][j] = GetGaussIntegral(i, j, Counts);
            int_matrix[j][i] = int_matrix[i][j];
        }
    }
    ShowIntMatrix(int_matrix);
    const str file_name = "params/SmartNormMap" + to_string(Counts) + ".txt";
    WriteIntMatrix(file_name, int_matrix);
    return int_matrix;
}

double DalitzMCIntegral::normal_pdf(const double x, const double m,
                                    const double s) {
    return inv_sqrt_2pi/s*exp(-0.5*pow((x-m)/s, 2));
}

double DalitzMCIntegral::GetGaussIntegral(const unsigned resnum,
                                          const uint64_t& nc) const {
    const uint64_t Counts = nc ? nc : m_ncounts;
    DStrip* strip = m_model->GetResAreas()[resnum];
    vectd mABsqV;
    vectd mACsqV;
    vectd mBCsqV;
    if (GetGaussPoints(&mABsqV, &mACsqV, &mBCsqV, Counts, strip)) return 0;
    const vectd& xV = strip->type == 1 ?
                mABsqV : (strip->type == 2 ? mACsqV : mBCsqV);
    double res = 0;
    const double mean = strip->mean();
    const double sigm = strip->sigma();
    for (uint64_t count=0; count < Counts; count++) {
        if (count && !(count%100000)) cout << count << " counts" << endl;
        res += norm(m_model->GetResVal(mABsqV[count], mACsqV[count], resnum)) /
                normal_pdf(xV[count], mean, sigm);
    }
    return res/Counts;
}

compld DalitzMCIntegral::GetGaussIntegral(const unsigned resnum1,
                                          const unsigned resnum2,
                                          const uint64_t& nc) const {
    const uint64_t Counts = nc ? nc : m_ncounts;
    DStrip* strip1 = m_model->GetResAreas()[resnum1];
    DStrip* strip2 = m_model->GetResAreas()[resnum2];
    vectd mABsqV;
    vectd mACsqV;
    vectd mBCsqV;
    compld res = 0;

    if (strip1->type == strip2->type) {
        const int rtype = strip1->type;
        const double mean = 0.5*(strip1->mean() + strip2->mean());
        const double sigm =
                std::max(0.5*std::fabs(strip1->mean()-strip2->mean()),
                         std::max(strip1->sigma(), strip2->sigma()));
        DStrip strip(mean-sigm, mean+sigm, rtype);
        if (GetGaussPoints(&mABsqV, &mACsqV, &mBCsqV, Counts, &strip)) return 0;
        const vectd& xV = rtype == 1 ? mABsqV : (rtype == 2 ? mACsqV : mBCsqV);
        for (uint64_t count=0; count < Counts; count++) {
            if (count && !(count%100000)) cout << count << " counts" << endl;
            const double gdenc = normal_pdf(xV[count], mean, sigm);
            const compld res1 = m_model->GetResVal(mABsqV[count],
                                                   mACsqV[count], resnum1);
            const compld res2 = m_model->GetResVal(mABsqV[count],
                                                   mACsqV[count], resnum2);
            res += res1*conj(res2)/gdenc;
        }
    } else {
        if (GetGaussPoints(&mABsqV, &mACsqV, &mBCsqV, Counts, strip1, strip2))
            return 0;
        const vectd& xV = strip1->type == 1 ?
                    mABsqV : (strip1->type == 2 ? mACsqV : mBCsqV);
        const vectd& yV = strip2->type == 1 ?
                    mABsqV : (strip2->type == 2 ? mACsqV : mBCsqV);
        const double meanX = strip1->mean();
        const double sigmX = strip1->sigma();
        const double meanY = strip2->mean();
        const double sigmY = strip2->sigma();
        for (uint64_t count=0; count < Counts; count++) {
            if (!(count%100000)) cout << count << " counts" << endl;
            const double gdenc = normal_pdf(xV[count], meanX, sigmX) *
                    normal_pdf(yV[count], meanY, sigmY);
            const compld res1 =
                    m_model->GetResVal(mABsqV[count], mACsqV[count], resnum1);
            const compld res2 =
                    m_model->GetResVal(mABsqV[count], mACsqV[count], resnum2);
            res += res1*conj(res2)/gdenc;
        }
    }
    return res*(1./Counts);
}

int DalitzMCIntegral::WriteIntMatrix(const str& fname,
                                     const vector<vectcd>& vals) const {
    ofstream ofile(fname.c_str(), ofstream::out);
    if (!ofile.is_open()) {
        cout << "Can't create file " << fname << endl;
        return -1;
    }
    const unsigned size = vals.size();
    ofile << "Size: " << size << endl;
    for (unsigned i=0; i < size; i++) {
        for (unsigned j=0; j < size; j++) {
            ofile << i << " " << j << " " << vals[i][j] << endl;
        }
    }
    ofile.close();
    return 0;
}

int DalitzMCIntegral::WriteIntMatrix(const str& fname,
                                     const vector<vectcd>& vals,
                                     const vector<vectcd>& errs) const {
    ofstream ofile(fname.c_str(), ofstream::out);
    if (!ofile.is_open()) {
        cout << "Can't create file " << fname << endl;
        return -1;
    }
    const unsigned size = vals.size();
    ofile << "Size: " << size << endl;
    for (unsigned i=0; i < size; i++) {
        for (unsigned j=0; j < size; j++) {
            ofile << i << " " << j << " " << vals[i][j]
                     << " +- " << errs[i][j] << endl;
        }
    }
    ofile.close();
    return 0;
}

void DalitzMCIntegral::ShowIntMatrix(const vector<vectcd>& vals) const {
    const unsigned size = vals.size();
    cout << "Size: " << size << endl;
    cout << std::setprecision(2) << std::scientific;
    for (unsigned i=0; i < size; i++) {
        for (unsigned j=0; j < size; j++) {
            cout << vals[i][j] << " ";
        }
        cout << endl;
    }
}

void DalitzMCIntegral::ShowIntMatrix(const vector<vectcd>& vals,
                                     const vector<vectcd>& errs) const {
    const unsigned size = vals.size();
    cout << "Size: " << size << endl;
    cout << std::setprecision(2) << std::scientific;
    for (unsigned i=0; i < size; i++) {
        for (unsigned j=0; j < size; j++) {
            cout << i << " " << j << " " << vals[i][j]
                    << " +- " << errs[i][j] << endl;
        }
    }
}
