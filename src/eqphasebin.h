/** Copyright 2017 Vitaly Vorobyev
 ** @file eqphasebin.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_EQPHASEBIN_H_
#define SRC_EQPHASEBIN_H_

#include "./absdalitzmodel.h"

///
/// \brief The EqPhaseBin class. Calculates equal phase bin number
/// for a given model.
///
class EqPhaseBin {
 public:
    EqPhaseBin(const AbsDalitzModel* model, const unsigned NBins);
    int Bin(const double& mABsq, const double& mACsq);
    double delta(const double& mABsq, const double& mACsq);

 private:
    double lEdge(const int i);
    double rEdge(const int i);

    const AbsDalitzModel* m_model;
    unsigned m_nbins;

    double m_del_min;
    double m_del_max;
};

#endif  // SRC_EQPHASEBIN_H_
