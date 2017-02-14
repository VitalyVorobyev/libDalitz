/** Copyright 2017 Vitaly Vorobyev
 ** @file bwwidth.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/bwwidth.h"

#include <cmath>

using std::pow;
using std::sqrt;

BWWidth::BWWidth(const double& G0, const double &m,
                 const double &p0, const int mom) :
    AbsVarWidth(G0, m, p0), m_mom(mom) {
    m_precalc = G0*m/pow(p0, 2*mom+1);
    m_ff = new BlattWeisskopf(mom, p0, BlattWeisskopf::FFResonance);
}

double BWWidth::operator()(const double& s, const double& p) const {
    const double ff = (*m_ff)(p);
    return m_precalc*ff*ff*pow(p, 2*m_mom+1)/sqrt(s);
}
