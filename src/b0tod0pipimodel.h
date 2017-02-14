/** Copyright 2017 Vitaly Vorobyev
 ** @file b0tod0pipimodel.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_B0TOD0PIPIMODEL_H_
#define SRC_B0TOD0PIPIMODEL_H_

#include "./symdalitzmodel.h"

/**
 * @brief The B0toD0pipiModel class. Implementation of B0 -> D0 pi+ pi- decay
 * amplitude. Belle and LHCb models are available.
/// The Belle one is from A. Kuzmin et al. (Belle Collaboration) Phys. Rev. D 76, 012006 (2007)
/// The LHCb  one is from R. Aaij   et al. (LHCb Collaboration)  Phys. Rev. D 92, 032002 (2015)
 */
class B0toD0pipiModel : public SymDalitzModel {
 public:
    explicit B0toD0pipiModel(const int type = LHCb);
    B0toD0pipiModel(const double& mB, const double& mD,
                    const double& mpi, const int type);

    static const int Belle;
    static const int LHCb;

 private:
    /**
     * @brief InitBelleModel. // ** A. Kuzmin et al. (Belle Collaboration)
     * Phys. Rev. D 76, 012006 – Published 30 July 2007
     */
    void InitBelleModel(void);
    /**
     * @brief InitLHCbModel. R. Aaij   et al. (LHCb Collaboration)
     * Phys. Rev. D 92, 032002 (2015)
     */
    void InitLHCbModel(void);
    int m_type;

    static const double m_B0_Mass;
    static const double m_D0_Mass;
    static const double m_PI_Mass;
    static const double dtr;
};

#endif  // SRC_B0TOD0PIPIMODEL_H_
