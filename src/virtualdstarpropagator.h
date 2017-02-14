/** Copyright 2017 Vitaly Vorobyev
 ** @file virtualdstarpropagator.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_VIRTUALDSTARPROPAGATOR_H_
#define SRC_VIRTUALDSTARPROPAGATOR_H_

#include <complex>

#include "./abspropagator.h"

///
/// \brief The VirtualDstarPropagator class
///
class VirtualDstarPropagator : public AbsPropagator {
 public:
    VirtualDstarPropagator(const double& beta1, const double& beta2);

    std::complex<double> operator()(const double& s,
                                    const double& p = 0) const;
 private:
    double m_b1;
    double m_b2;
};

#endif  // SRC_VIRTUALDSTARPROPAGATOR_H_
