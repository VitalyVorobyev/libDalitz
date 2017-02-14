/** Copyright 2017 Vitaly Vorobyev
 ** @file gswidth.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/gswidth.h"
#include <cmath>

GSWidth::GSWidth(const double &G0, const double &m, const double &p0) :
    AbsVarWidth(G0, m, p0), m_precalc(G0*m/std::pow(p0, 3)) {}

double GSWidth::operator()(const double& s, const double& p) const {
    return m_precalc*std::pow(p, 3)/std::sqrt(s);
}
