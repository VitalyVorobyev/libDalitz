/** Copyright 2017 Vitaly Vorobyev
 ** @file smalldpbin.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
  **/

#include "../src/smalldpbin.h"

#include <cmath>

typedef std::complex<double> compld;

SmallDPBin::SmallDPBin(const double &mAB, const double &mAC) :
    m_val(0), m_valb(0),
    m_wght(0), m_dph(0), m_p(0), m_pb(0),
    m_h(0), m_gsize(0), m_norm(1),
    m_dpsize(0), m_dpmin(0),
    m_mAB(mAB), m_mAC(mAC),
    m_ABbin(GetBin(m_mAB)),
    m_ACbin(GetBin(m_mAC)) {}

SmallDPBin::SmallDPBin(AbsDalitzModel* model, const double& mAB,
                       const double& mAC) :
    SmallDPBin(mAB, mAC) {
    SetModel(model);
    Calc();
}

SmallDPBin::SmallDPBin(const SmallDPBin& x) {
    *this = x;
}

void SmallDPBin::SetModel(AbsDalitzModel* model) {
    m_model  = model;
    m_dpmin  = m_model->mABsq_min();
    m_dpsize = m_model->mABsq_max() - m_dpmin;
}

void SmallDPBin::SetGridSize(const int gsize) {
    m_gsize = gsize;
    m_h = m_dpsize / m_gsize;
}

void SmallDPBin::SetGridStep(const double& gstep) {
    m_h = gstep;
    m_gsize = static_cast<int>(m_dpsize / m_h);
}

void SmallDPBin::SetPoint(const double& mAB, const double& mAC) {
    if (!m_model->IsInPlot(mAB, mAC)) return;
    m_mAB = mAB;
    m_mAC = mAC;
    m_ABbin = GetBin(m_mAB);
    m_ACbin = GetBin(m_mAC);
    Calc();
}

void SmallDPBin::SetBin(const int i, const int j) {
    m_ABbin = i;
    m_ACbin = j;
    m_mAB = GetVal(m_ABbin);
    m_mAC = GetVal(m_ACbin);
    if (!m_model->IsInPlot(m_mAB, m_mAC)) return;
    Calc();
}

void SmallDPBin::SetNorm(const double& norm) {
    m_norm = norm;
}

bool SmallDPBin::IsInBin(const double& mAB, const double& mAC) const {
    if (std::fabs(mAB-m_mAB) > m_h) return false;
    if (std::fabs(mAC-m_mAC) > m_h) return false;
    return true;
}

double SmallDPBin::Phase(void) const {return arg(m_val);}
double SmallDPBin::Amp(void) const {return sqrt(m_p);}
double SmallDPBin::Ampb(void) const {return sqrt(m_pb);}
double SmallDPBin::P(void) const {return m_p;}
double SmallDPBin::Pb(void) const {return m_pb;}
compld SmallDPBin::CAmp(void) const {return m_val;}
compld SmallDPBin::CAmpb(void) const {return m_valb;}
int SmallDPBin::GSize(void) const {return m_gsize;}
double SmallDPBin::Weight(void) const {return m_wght;}

void SmallDPBin::GetCurrentPoint(double* mAB, double* mAC) const {
    *mAB = m_mAB;
    *mAC = m_mAC;
}

int SmallDPBin::GetBin(const double& x) const {
    if (x <= 0 || x > (m_dpmin + m_dpsize) || x < m_dpmin) return -1;
    return (x - m_dpmin) / m_h;
}

double SmallDPBin::GetVal(const int bin) const {
    if (bin < 0 || bin >= m_gsize) return 0;
    return m_dpmin + m_h * (static_cast<double>(bin) + 0.5);
}

int SmallDPBin::Calc(void) {
    if (!m_model->IsInPlot(m_mAB, m_mAC)) return -1;
    m_val = m_model->Amp(m_mAB, m_mAC);
    m_valb = m_model->Amp(m_mAC, m_mAB);
    m_dph = std::arg(m_valb / m_val);
    m_p = std::norm(m_val);
    m_pb = std::norm(m_valb);
    m_wght = std::sqrt(m_p * m_pb) * (m_p + m_pb);
    return 0;
}
