/** Copyright 2017 Vitaly Vorobyev
 ** @file gounarissakurai.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/gounarissakurai.h"

#include <cmath>
#include <iostream>

#include "./constwidth.h"
#include "./gswidth.h"

using std::log;

const double GounarisSakurai::mpi = 0.13957018;

GounarisSakurai::GounarisSakurai(const double &G0, const double &m,
                                 const double &p0, const bool constwidth) :
    AbsPropagator(m, p0), m_const_width(constwidth),
    m_width(m_const_width ?
                reinterpret_cast<AbsVarWidth*>(new ConstWidth(G0)) :
                reinterpret_cast<AbsVarWidth*>(new GSWidth(G0, m, p0))) {
    g();
}

std::complex<double> GounarisSakurai::operator()(const double& s,
                                                 const double& p) const {
    const double& mr  = m();
    const double mrsq = mr*mr;
    const double& G0  = m_width->G0();
    const double Ggs  = (*m_width)(s, p);
    const auto imone = std::complex<double>(0, 1);
    return mrsq*(1.+G0*m_g/mr)/(mrsq-s+f(s, p)-imone*mr*Ggs);
}

GounarisSakurai::~GounarisSakurai() {
    m_const_width ? delete reinterpret_cast<ConstWidth*>(m_width) :
                    delete reinterpret_cast<GSWidth*>(m_width);
}

double GounarisSakurai::g(void) {
    const double  mpisq = mpi*mpi;
    const double& pi    = M_PI;
    const double  p0sq  = p0()*p0();
    const double& mr    = m();
    return m_g = 3./pi*mpisq/(p0sq)*log((mr+2.*p0())/2.*mpi) + mr/(2.*pi*p0())
            - mpisq*mr/(pi*p0sq*p0());
}

int GounarisSakurai::hhder(double* h, double* hder, const double& p) const {
    const double& mr   = m();
    const double coeff = 2./M_PI;

    const double var = (mr+2.*p)/(2.*mpi);
    if (var <= 0) return -1;

    const double logi = log(var);
    *h    = coeff*p/mr*logi;
    *hder = coeff*mpi*p/(mr*mr)*(1./(mr+2.*p)-logi/mr);
    return 0;
}

double GounarisSakurai::h(const double& s, const double& p) const {
    const double sqrts = sqrt(s);
    const double coeff = 2./M_PI;

    const double var = (sqrts+2.*p)/(2.*mpi);
    if (var <= 0) return -1;

    const double logi = log(var);
    return coeff*p/sqrts*logi;
}

double GounarisSakurai::f(const double& s, const double& p) const {
    double hr, hrder;
    hhder(&hr, &hrder, p);
    const double hs = h(s, p);

    const double& G0  = m_width->G0();
    const double mrsq = m()*m();
    const double p0sq = p0()*p0();

    return G0*mrsq/(p0()*p0sq)*((hs-hr)*p*p + (mrsq-s)*p0sq*hrder);
}
