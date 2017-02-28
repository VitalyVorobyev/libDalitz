/** Copyright 2017 Vitaly Vorobyev
 ** @file modelintegral.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/modelintegral.h"

#include <sstream>
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

    double mmMax, mmMin;
    const double M_min = m_model->mABsq_min();
    const double M_max = m_model->mABsq_max();
    const double dm = (M_max - M_min) / m_gsize;
    const str file_name_binning = "params/" + label + "_binning.txt";
    std::ofstream file_binning(file_name_binning);
    file_binning << "GridSize " << m_gsize << endl;
    file_binning << m_model->AsText() << endl;
    double mp = M_min;
    for (unsigned i=0; i < m_gsize; i++) {
        mp += dm;
        m_model->mABsqRange_AC(mp, &mmMin, &mmMax);
        if (mp > mmMax) break;
        for (double mm = mp; mm < mmMax; mm += dm) {
            if (!m_model->IsInPlot(mp, mm)) continue;
            int bin = abs(eqbin.Bin(mp, mm));
            file_binning << mp << " " << mm << " "
                         << m_model->m3sq(mp, mm) << " " << bin << endl;
            double delta, P, Pbar;
            if (PPbarDelta(mp, mm, &P, &Pbar, &delta)) continue;
            bin -= 1;

            C->at(bin)  += sqrt(P*Pbar)*std::cos(delta);
            S->at(bin)  += sqrt(P*Pbar)*std::sin(delta);
            Kp->at(bin) += P;
            Kn->at(bin) += Pbar;
        }
    }
    file_binning.close();

    double norm = 0;
    for (unsigned i=0; i < m_nbins; i++) {
        const double rmp = sqrt(Kp->at(i)*Kn->at(i));
        C->at(i) /= rmp;
        S->at(i) /= rmp;
        norm += Kp->at(i) + Kn->at(i);
    }

    cout << label << endl;
    const str file_name_binned_params("params/test.txt");
    cout << file_name_binned_params << endl;
    std::ofstream file_binned_params(file_name_binned_params);
    cout << "Norm     = " << norm << endl;
    unsigned i = 0;
    std::stringstream out;
    for (auto KpIt = Kp->begin(), KnIt = Kn->begin(),
              CIt = C->begin(), SIt = S->begin(); i < m_nbins;
         i++, KpIt++, KnIt++, CIt++, SIt++) {
        *KpIt /= norm; *KnIt /= norm;

        out.str("");
        out << i+1 << ": C = "  << *CIt << ", S = "  << *SIt
            << ", Kp = " << *KpIt << ", Kn = " << *KnIt
            << ", Q = "  << pow(*CIt, 2)+pow(*SIt, 2) << endl;
        cout << out.str();
        file_binned_params << out.str();
    }
    file_binned_params.close();
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
