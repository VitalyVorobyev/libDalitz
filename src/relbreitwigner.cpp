/** Copyright 2017 Vitaly Vorobyev
 ** @file relbreitwigner.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/relbreitwigner.h"

#include <iostream>

#include "./bwwidth.h"
#include "./constwidth.h"
#include "./flattewidth.h"

RelBreitWigner::RelBreitWigner(const double &G0, const double &m,
                               const double &p0, const int mom,
                               const int wtype) :
    AbsPropagator(m, p0), m_wtype(wtype) {
    switch (m_wtype) {
    case VarWType::Const:
        m_width = new ConstWidth(G0);
        break;
    case VarWType::BW:
        m_width = new BWWidth(G0, m, p0, mom);
        break;
    case VarWType::Flatte:
        m_width = new FlatteWidth(m);
        break;
    default:
        std::cout << "RelBreitWigner: wrong VarWType " << m_wtype << std::endl;
        break;
    }
}

std::complex<double> RelBreitWigner::operator()(const double& s,
                                                const double& p) const {
    const std::complex<double> ione(0, 1);
    const double& mass = m();
    return 1. / (mass * mass - s - ione * mass * (*m_width)(s, p));
}

RelBreitWigner::~RelBreitWigner() {
    delete m_width;
}
