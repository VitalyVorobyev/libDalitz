/** Copyright 2017 Vitaly Vorobyev
 ** @file virtualresff.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "../src/virtualresff.h"
#include <cmath>

VirtualResFF::VirtualResFF(const double& _r, const double& _p0) :
  FormFactor(_r, _p0) {}

double VirtualResFF::operator()(const double& p) const {
    return std::exp(-1. * std::fabs(p - p0()));
}
