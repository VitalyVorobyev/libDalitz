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

#include "./abssymdalitzmodel.h"
#include "./dalitzmodel.h"

/**
 * @brief The B0toD0pipiModel class. Implementation of B0 -> D0 pi+ pi- decay
 * amplitude. Belle and LHCb models are available.
 * The Belle one is from A. Kuzmin et al. (Belle Collaboration)
 * Phys. Rev. D 76, 012006 (2007)
 * The LHCb  one is from R. Aaij   et al. (LHCb Collaboration)
 * Phys. Rev. D 92, 032002 (2015)
 */
class B0toD0pipiModel : public DalitzModel, public AbsSymDalitzModel {
 public:
    explicit B0toD0pipiModel(int type = LHCb);
    B0toD0pipiModel(double mB, double mD, double mpi, int type);

    static constexpr int Belle = 0;
    static constexpr int LHCb = 1;

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

    static constexpr double m_B0_Mass = 5.279;
    static constexpr double m_D0_Mass = 1.865;
    static constexpr double m_PI_Mass = 0.139568;
    static constexpr double dtr = M_PI / 180.;
};

#endif  // SRC_B0TOD0PIPIMODEL_H_
