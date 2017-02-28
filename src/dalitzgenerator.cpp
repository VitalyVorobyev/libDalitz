/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzgenerator.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/dalitzgenerator.h"

#include <iostream>
#include <fstream>

using std::ofstream;

typedef DalitzGenerator DGen;
typedef std::uniform_real_distribution<double> unifdist;
typedef std::vector<double> vectd;

using std::cout;
using std::cerr;
using std::endl;

DalitzGenerator::DalitzGenerator(const AbsDalitzModel *model) :
    RandomDalitzPoint(model), m_ntries(1e8), m_maj_counts(1e4),
    m_maj(-1), m_model(model) {
    SetMajorant(CalcMajorant());
}

void DGen::SetMajorant(const double& p) {
    m_maj = p;
    cout << "DalitzGenerator: majorant is set to " << m_maj << endl;
}

double DGen::CalcMajorant(void) const {
    cout << "CalcMajorant..." << endl;
    double cur_val;
    double maj = 0;
    double mABsq, mACsq;
    for (unsigned i=0; i < m_maj_counts; i++) {
        GetPoint(&mABsq, &mACsq);
        cur_val = m_model->P(mABsq, mACsq);
        if (cur_val > maj) maj = cur_val;
    }
    maj *= 1.05;
    return maj;
}

int DGen::Generate(double *mABsq, double *mACsq, unifdist* dist) const {
    uint64_t tries = 0;
    if (dist == nullptr) *dist = unifdist(0., m_maj);
    while (tries++ < m_ntries) {
        GetPoint(mABsq, mACsq);
        if ((*dist)(*re) < m_model->P(*mABsq, *mACsq)) break;
    }
    if (tries == m_ntries) {
        cerr << "DalitzGenerator::Generate: tries limit exceeded!" << endl;
        return 0;
    }
    return tries;
}

int DGen::Generate(const uint64_t NEv, vectd* mABv, vectd* mACv,
                   const bool silent) const {
    mABv->resize(NEv); mACv->resize(NEv);
    uint64_t tries = 0;
    uint64_t Evtn = 0;
    unifdist unifMaj(0., m_maj);
    cout << "Generating " << NEv << " events..." << endl;
    for (auto itAB = mABv->begin(), itAC = mACv->begin();
         itAB != mABv->end() && itAC != mACv->end(); itAB++, itAC++) {
        const unsigned ntries = Generate(&*itAB, &*itAC, &unifMaj);
        if (ntries == 0) break;
        tries += ntries;
        Evtn++;
    }
    if (!silent) {
        const double Eff = static_cast<double>(Evtn)/tries;
        if (tries == m_ntries) {
            cout << "DalitzGenerator: tries limit exceeded!";
            cout << "  efficiency: " << Eff << endl;
            cout << "  events generated: " << Evtn << endl;
            cout << "  tries limit: " << m_ntries << endl;
            return -1;
        }
        cout << "Done!. Efficiency " << Eff << endl;
    }
    return 0;
}

int DGen::WriteDDist(const uint64_t& NEv, const std::string& fname) const {
    ofstream file(fname, ofstream::out);
    if (!file.is_open()) {
        cerr << "Can't open file " << fname << endl;
        return -1;
    }
    file << NEv << " events" << endl << m_model->AsText();
    uint64_t eve_counter = 0;
    double mAB, mAC;
    unifdist unifMaj(0., m_maj);
    while (eve_counter++ < NEv) {
        if (Generate(&mAB, &mAC, &unifMaj) == 0) break;
        file << mAB << " " << mAC << endl;
    }
    file.close();
    return eve_counter != NEv;
}
