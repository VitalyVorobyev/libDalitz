/** Copyright 2017 Vitaly Vorobyev
 ** @file gswidth.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_GSWIDTH_H_
#define SRC_GSWIDTH_H_

#include "./absvarwidth.h"

class GSWidth : public AbsVarWidth {
 public:
    GSWidth(const double& G0, const double& m, const double& p0);

    double operator()(const double& s, const double& p) const;

 private:
    double m_precalc;
};

#endif  // SRC_GSWIDTH_H_
