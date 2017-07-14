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
#include <algorithm>

typedef std::vector<double> vectd;
typedef std::string str;

using std::cout;
using std::endl;
using std::abs;

ModelIntegral::ModelIntegral(const AbsSymDalitzModel *model,
                             const unsigned gridsize)  :
    m_model(model), m_gsize(gridsize) {}

int ModelIntegral::Calculate(const str& label,
                             vectd* C, vectd* S, vectd* Kp, vectd* Kn) {
    const unsigned nbins = m_model->nbins();
    C->resize(nbins, 0);
    S->resize(nbins, 0);
    Kp->resize(nbins, 0);
    Kn->resize(nbins, 0);

    const double M_min = m_model->mABsq_min();
    const double M_max = m_model->mABsq_max();
    const double dm = (M_max - M_min) / m_gsize;
//    const str file_name_binning = "params/" + label + "_binning.txt";
//    std::ofstream file_binning(file_name_binning);
//    file_binning << "GridSize " << m_gsize << endl;
//    file_binning << m_model->AsText() << endl;
    double mp = M_min;
    for (unsigned i=0; i < m_gsize; i++) {
        mp += dm;
        double mmMax, mmMin;
        m_model->mABsqRange_AC(mp, &mmMin, &mmMax);
        if (mp > mmMax) break;
        for (double mm = std::max(mp, mmMin); mm < mmMax; mm += dm) {
            if (!m_model->IsInPlot(mp, mm)) continue;
            double delta, p, pb;
            m_model->ppdelt(mp, mm, &p, &pb, &delta);
            int bidx = m_model->bin(delta, false) - 1;
//            file_binning << mp << " " << mm << " "
//                         << m_model->m3sq(mp, mm) << " " << bidx+1 << endl;
            const double sqppb = sqrt(p*pb);
            C->at(bidx)  += sqppb*std::cos(delta);
            S->at(bidx)  += sqppb*std::sin(delta);
            Kp->at(bidx) += p;
            Kn->at(bidx) += pb;
        }
    }
//    file_binning.close();

    double norm = 0;
    for (unsigned i=0; i < nbins; i++) {
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
              CIt = C->begin(), SIt = S->begin(); i < nbins;
         i++, KpIt++, KnIt++, CIt++, SIt++) {
        *KpIt /= norm; *KnIt /= norm;

        out.str("");
        out << i+1 << ": C = "  << *CIt << ", S = "  << *SIt
            << ", Kp = " << *KpIt << ", Kn = " << *KnIt
//            << ", Q = "  << pow(*CIt, 2)+pow(*SIt, 2)
            << endl;
        cout << out.str();
        file_binned_params << out.str();
    }
    file_binned_params.close();
    return 0;
}
