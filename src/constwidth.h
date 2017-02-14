/** Copyright 2017 Vitaly Vorobyev
 ** @file constwidth.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_CONSTWIDTH_H_
#define SRC_CONSTWIDTH_H_

#include "./absvarwidth.h"

class ConstWidth : public AbsVarWidth {
 public:
    explicit ConstWidth(const double& G0);

    double operator()(const double& s = 0,
                      const double& p = 0) const {return G0();}
};

#endif  // SRC_CONSTWIDTH_H_
