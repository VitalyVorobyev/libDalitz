/** Copyright 2017 Vitaly Vorobyev
 ** @file modelintegral.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_MODELINTEGRAL_H_
#define SRC_MODELINTEGRAL_H_

#include <string>
#include <vector>

#include "./absdalitzmodel.h"

///
/// \brief The ModelIntegral class. Class for calculation binned Dalitz plot
/// parameters using a gived decay amplitude model.
///
class ModelIntegral {
 public:
    ModelIntegral(const AbsDalitzModel* model, const unsigned NBins = 8,
                  const unsigned gridsize = 1000);

    void SetGridSize(const unsigned gsize) {m_gsize = gsize;}
    void SetNBins(const unsigned nbins) {m_nbins = nbins;}
    int Calculate(const std::string& label,
                  std::vector<double>* C, std::vector<double>* S,
                  std::vector<double>* K, std::vector<double>* Kb);

 private:
    int PPbarDelta(const double& mABsq, const double& mACsq,
                   double* P, double* Pbar, double* delta);
    const AbsDalitzModel* m_model;
    unsigned m_nbins;
    unsigned m_gsize;
};

#endif  // SRC_MODELINTEGRAL_H_
