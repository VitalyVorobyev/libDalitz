/** Copyright 2017 Vitaly Vorobyev
 ** @file kspipimodel.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_KSPIPIMODEL_H_
#define SRC_KSPIPIMODEL_H_

#include <cmath>

#include "./abssymdalitzmodel.h"
#include "./dalitzmodel.h"

///
/// \brief Implementation of D0 -> Ks0 pi+ pi- Belle 2010 decay model.
/// Parameters are taken from A. Poluektov et al.
/// Phys. Rev. D 81, 112002 – Published 16 June 2010.
/// See also S. Kopp et al. (CLEO Collaboration)
/// Phys. Rev. D 63, 092001 – Published 9 April 2001
///
class KspipiModel :
        public DalitzModel,
        public AbsSymDalitzModel {

    static constexpr double m_D0_Mass = 1.865;
    static constexpr double m_Ks0_Mass = 0.497611;
    static constexpr double m_PI_Mass = 0.13957018;
    static constexpr double dtr = M_PI / 180.;

 public:
   KspipiModel(void);
   KspipiModel(double md, double mks, double mpi);
};

#endif  // SRC_KSPIPIMODEL_H_
