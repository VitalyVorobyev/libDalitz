/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzvetodst2010.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#include "../src/dalitzvetodst2010.h"
#include <cmath>

DalitzVetoDst2010::DalitzVetoDst2010(const double cut) :
    m_mdst(2.0103), m_cut(cut) {}

bool DalitzVetoDst2010::operator()(const double& AB, const double& BC) const {
    return std::fabs(std::sqrt(AB)-m_mdst) < m_cut;
}
