/** Copyright 2017 Vitaly Vorobyev
 ** @file virtualdstarpropagator2.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_VIRTUALDSTARPROPAGATOR2_H_
#define SRC_VIRTUALDSTARPROPAGATOR2_H_

#include <complex>

#include "./abspropagator.h"

class VirtualDstarPropagator2 : public AbsPropagator {
 public:
    VirtualDstarPropagator2();
    std::complex<double> operator()(const double& s, const double& p) const;

 private:
    double m_r;
};

#endif  // SRC_VIRTUALDSTARPROPAGATOR2_H_
