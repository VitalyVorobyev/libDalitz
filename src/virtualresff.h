/** Copyright 2017 Vitaly Vorobyev
 ** @file virtualresff.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_VIRTUALRESFF_H_
#define SRC_VIRTUALRESFF_H_

#include "./formfactor.h"

///
/// \brief The VirtualResFF class
///
class VirtualResFF : public FormFactor {
 public:
    VirtualResFF(const double& _r, const double& _p0);
    double operator()(const double& p) const;
};

#endif  // SRC_VIRTUALRESFF_H_
