/** Copyright 2017 Vitaly Vorobyev
 ** @file buggproparator.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_BUGGPROPAGATOR_H_
#define SRC_BUGGPROPAGATOR_H_

#include <complex>

#include "./abspropagator.h"
#include "./buggwidth.h"

class BuggPropagator : public AbsPropagator {
 public:
    BuggPropagator(void);
    ~BuggPropagator(void);

    std::complex<double> operator()(const double& s,
                                    const double& p = 0) const;

 private:
    BuggWidth* m_width;
};

#endif  // SRC_BUGGPROPAGATOR_H_
