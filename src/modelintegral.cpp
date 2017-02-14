/** Copyright 2017 Vitaly Vorobyev
 ** @file modelintegral.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/modelintegral.h"

#include <fstream>
#include <cmath>
#include <iostream>

#include "./eqphasebin.h"

typedef std::vector<double> vectd;
typedef std::string str;

using std::cout;
using std::endl;
using std::abs;

ModelIntegral::ModelIntegral(const AbsDalitzModel *model,
                             const unsigned NBins,
                             const unsigned gridsize)  :
    m_model(model), m_nbins(NBins), m_gsize(gridsize) {}

int ModelIntegral::Calculate(const str& label,
                             vectd* C, vectd* S, vectd* Kp, vectd* Kn) {
    C->resize(m_nbins, 0);
    S->resize(m_nbins, 0);
    Kp->resize(m_nbins, 0);
    Kn->resize(m_nbins, 0);
    EqPhaseBin eqbin(m_model, m_nbins);

    double mmMax;
    double mmMin;
    const double M_min = m_model->mABsq_min();
    const double M_max = m_model->mABsq_max();
    const double dm = (M_max - M_min)/m_gsize;
    const str fname = "params/" + label + "_binning.txt";
    std::ofstream ofile(fname.c_str());
    ofile << "GridSize " << m_gsize << endl;
    ofile << m_model->mM() << " -> " << m_model->mA();
    ofile << " " << m_model->mB();
    ofile << " " << m_model->mC() << endl;
    ofile << "mAB: " << m_model->mABsq_min() << " "
          << m_model->mABsq_max() << endl;
    ofile << "mAC: " << m_model->mACsq_min() << " "
          << m_model->mACsq_max() << endl;
    ofile << "mBC: " << m_model->mBCsq_min() << " "
          << m_model->mBCsq_max() << endl;
    ofile << m_model->ABaxis() << endl;
    ofile << m_model->ACaxis() << endl;
    ofile << m_model->BCaxis() << endl;
    double mp = M_min;
    for (unsigned i=0; i < m_gsize; i++) {
        mp += dm;
        m_model->mABsqRange_AC(mp, &mmMin, &mmMax);
        if (mp > mmMax) break;
        double mm = mp;
        for (int j=0; mm < mmMax; j++) {
            mm += dm;
            if (!m_model->IsInPlot(mp, mm)) continue;
            int bin = abs(eqbin.Bin(mp, mm));
            ofile << mp << " " << mm << " " << m_model->GetmBCsq(mp, mm)
                  << " " << bin << endl;
            double delta, P, Pbar;
            if (PPbarDelta(mp, mm, &P, &Pbar, &delta)) continue;
            bin -= 1;

            C->at(bin)  += sqrt(P*Pbar)*std::cos(delta);
            S->at(bin)  += sqrt(P*Pbar)*std::sin(delta);
            Kp->at(bin) += P;
            Kn->at(bin) += Pbar;
        }
    }
    ofile.close();

    double norm = 0;
    for (unsigned i=0; i < 8; i++) {
        const double rmp = sqrt(Kp->at(i)*Kn->at(i));
        C->at(i) /= rmp;
        S->at(i) /= rmp;
        norm += Kp->at(i) + Kn->at(i);
    }

    cout << "Hello, God! " << label << endl;
    const str fname1("params/test.txt");
    cout << fname1 << endl;
    std::ofstream ofile1(fname1.c_str());
    cout << "Norm     = " << norm << endl;
    for (unsigned i=0; i < m_nbins; i++) {
        Kp->at(i) /= norm;
        Kn->at(i) /= norm;

        cout << i+1;
        cout << ": C = "  << C->at(i);
        cout << ", S = "  << S->at(i);
        cout << ", Kp = " << Kp->at(i);
        cout << ", Kn = " << Kn->at(i);
        cout << ", Q = "  << C->at(i)*C->at(i)+S->at(i)*S->at(i) << endl;

        ofile1 << i+1;
        ofile1 << ": C = "  << C->at(i);
        ofile1 << ", S = "  << S->at(i);
        ofile1 << ", Kp = " << Kp->at(i);
        ofile1 << ", Kn = " << Kn->at(i);
        ofile1 << ", Q = "  << C->at(i)*C->at(i)+S->at(i)*S->at(i) << endl;
    }
    ofile1.close();
    return 0;
}

int ModelIntegral::PPbarDelta(const double& mABsq, const double& mACsq,
                              double* P, double* Pbar, double* delta) {
    std::complex<double> A  = m_model->Amp(mABsq, mACsq);
    std::complex<double> Ab = m_model->Amp(mACsq, mABsq);
    *delta = -(std::arg(Ab) - std::arg(A));
    if (std::isnan(*delta)) { cout << "delta: " << delta << endl; return -1;}
    *P = norm(A); *Pbar = norm(Ab);
    return 0;
}
