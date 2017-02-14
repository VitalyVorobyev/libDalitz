/** Copyright 2017 Vitaly Vorobyev
 ** @file symdalitzmodel.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_SYMDALITZMODEL_H_
#define SRC_SYMDALITZMODEL_H_

#include "./dalitzmodel.h"

class SymDalitzModel : public DalitzModel {
 public:
    SymDalitzModel(const double& mmo, const double& mcha, const double& mchb,
                   const double& delmin, const double& delmax);

    double delta(const double& mp, const double& mm);
    void PPbarDelta(const double& mp, const double& mm, double* P,
                    double* Pbar, double* delta);
    int GetBin(const double& mp, const double& mm);

    void SetNBins(const int nb) {nbins = nb; return;}
    int GetNBins(void) const {return nbins;}
 private:
    double del_min, del_max;
    int nbins;
};

#endif  // SRC_SYMDALITZMODEL_H_
