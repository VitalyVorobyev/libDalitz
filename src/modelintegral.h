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

#include "./abssymdalitzmodel.h"

/**
 * @brief The ModelIntegral class. Class for calculation binned Dalitz plot
 * parameters using a gived decay amplitude model.
 */
class ModelIntegral {
 public:
    ModelIntegral(const AbsSymDalitzModel* model,
                  const unsigned gridsize = 1000);

    void SetGridSize(const unsigned gsize) {m_gsize = gsize;}
    int Calculate(const std::string& label,
                  std::vector<double>* C, std::vector<double>* S,
                  std::vector<double>* K, std::vector<double>* Kb);

 private:
    const AbsSymDalitzModel* m_model;
    unsigned m_gsize;
};

#endif  // SRC_MODELINTEGRAL_H_
