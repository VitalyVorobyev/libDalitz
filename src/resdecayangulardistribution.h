/** Copyright 2017 Vitaly Vorobyev
 ** @file resdecayangulardistribution.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_RESDECAYANGULARDISTRIBUTION_H_
#define SRC_RESDECAYANGULARDISTRIBUTION_H_

///
/// \brief Angular pdf for intermidiate resonance decay in a three-body decay.
/// Final state particles are assumed to be scalars.
/// Resonance may be scalar, vector or tensor.
class ResDecayAngularDistribution {
 public:
    ResDecayAngularDistribution(const int spin, const double& mmo,
                                const double& mca, const double& mcb,
                                const double& mcc, const double& mres);
    double operator()(const double& mACsq, const double& mBCsq) const;

    void Set_mR(const double& mr) {m_mre_sq = mr * mr;}

    double mMotSq(void) const {return m_mmo_sq;}
    double mChASq(void) const {return m_mca_sq;}
    double mChBSq(void) const {return m_mcb_sq;}
    double mChCSq(void) const {return m_mcc_sq;}
    double mResSq(void) const {return m_mre_sq;}

    double mMot(void) const {return m_mmo;}
    double mChA(void) const {return m_mca;}
    double mChB(void) const {return m_mcb;}
    double mChC(void) const {return m_mcc;}

    static bool m_use_mRsq;

 private:
    int m_spin;

    double m_mmo;
    double m_mca;
    double m_mcb;
    double m_mcc;

    double m_mmo_sq;
    double m_mca_sq;
    double m_mcb_sq;
    double m_mcc_sq;
    double m_mre_sq;

    double m_pc_var1;
    double m_pc_var2;
    double m_pc_var3;
    double m_pc_var4;
    double m_pc_var5;

    void make_precalc_one(void);
    void make_precalc_all(void);
};

#endif  // SRC_RESDECAYANGULARDISTRIBUTION_H_
