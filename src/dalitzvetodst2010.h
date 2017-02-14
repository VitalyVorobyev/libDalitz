/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzvetodst2010.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_DALITZVETODST2010_H_
#define SRC_DALITZVETODST2010_H_

#include "./dalitzveto.h"

class DalitzVetoDst2010 : public DalitzVeto {
 public:
    explicit DalitzVetoDst2010(const double cut);
    bool operator()(const double& AB, const double& BC = 0) const;

 private:
    double m_mdst;
    double m_cut;
};

#endif  // SRC_DALITZVETODST2010_H_
