/** Copyright 2017 Vitaly Vorobyev
 ** @file nrpropagator.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_NRPROPAGATOR_H_
#define SRC_NRPROPAGATOR_H_

#include <complex>

#include "./abspropagator.h"

class NRPropagator : public AbsPropagator {
 public:
    explicit NRPropagator(const double& alpha);
    std::complex<double> operator()(const double& s,
                                    const double& p = 0) const;

 private:
    double m_alpha;
};

#endif  // SRC_NRPROPAGATOR_H_
