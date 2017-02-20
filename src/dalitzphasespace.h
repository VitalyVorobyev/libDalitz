/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzphasespace.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_DALITZPHASESPACE_H_
#define SRC_DALITZPHASESPACE_H_

#include <vector>
#include <complex>

#include "mylibs/libLinAl/lvect.h"

/**
 * @brief The DalitzPhaseSpace class. Comprehensive description of
 * a three-body decay kinematics
 */
class DalitzPhaseSpace {
 public:
    /**
     * @brief DalitzPhaseSpace. M -> ABC decay kinematics.
     * @param mM. Mass of the mother particle
     * @param mA. Mass of the daughter particle A.
     * @param mB. Mass of the daughter particle B.
     * @param mC. Mass of the daughter particle C.
     */
    DalitzPhaseSpace(const double& mM, const double& mA,
                     const double& mB, const double& mC);
    /**
     * @brief IsInPlot. Chack
     * @param mABsq
     * @param mACsq
     * @return true if the point is within allowed space, false otherwise
     */
    bool IsInPlot(const double& mABsq, const double& mACsq) const;
    ///
    /// \brief GetmBCsq
    /// \param mABsq
    /// \param mACsq
    /// \return
    ///
    double GetmBCsq(const double& mABsq, const double& mACsq) const;
    ///
    /// \brief mA
    /// \return
    ///
    double mA(void) const {return m_mA;}
    ///
    /// \brief mB
    /// \return
    ///
    double mB(void) const {return m_mB;}
    ///
    /// \brief mC
    /// \return
    ///
    double mC(void) const {return m_mC;}
    ///
    /// \brief mM
    /// \return
    ///
    double mM(void) const {return m_mM;}
    ///
    /// \brief mABsq_min
    /// \return
    ///
    double mABsq_min(void) const {return m_mABsq_min;}
    ///
    /// \brief mABsq_max
    /// \return
    ///
    double mABsq_max(void) const {return m_mABsq_max;}
    ///
    /// \brief mACsq_min
    /// \return
    ///
    double mACsq_min(void) const {return m_mACsq_min;}
    ///
    /// \brief mACsq_max
    /// \return
    ///
    double mACsq_max(void) const {return m_mACsq_max;}
    ///
    /// \brief mBCsq_min
    /// \return
    ///
    double mBCsq_min(void) const {return m_mBCsq_min;}
    ///
    /// \brief mBCsq_max
    /// \return
    ///
    double mBCsq_max(void) const {return m_mBCsq_max;}

    // * Ranges for squared inv masses * //
    /**
     * @brief mABsqRange_AC
     * @param mACsq
     * @param mABsqMin
     * @param mABsqMax
     * @return
     */
    int mABsqRange_AC(const double& mACsq, double* mABsqMin,
                      double* mABsqMax) const;
    /**
     * @brief mBCsqRange_AC
     * @param mACsq
     * @param mBCsqMin
     * @param mBCsqMax
     * @return
     */
    int mBCsqRange_AC(const double& mACsq, double* mBCsqMin,
                      double* mBCsqMax) const;
    /**
     * @brief mBCsqRange_AB
     * @param mABsq
     * @param mBCsqMin
     * @param mBCsqMax
     * @return
     */
    int mBCsqRange_AB(const double& mABsq,
                      double* mBCsqMin, double* mBCsqMax) const;
    /**
     * @brief mACsqRange_AB
     * @param mABsq
     * @param mACsqMin
     * @param mACsqMax
     * @return
     */
    int mACsqRange_AB(const double& mABsq,
                      double* mACsqMin, double* mACsqMax) const;
    /**
     * @brief mACsqRange_BC
     * @param mBCsq
     * @param mACsqMin
     * @param mACsqMax
     * @return
     */
    int mACsqRange_BC(const double& mBCsq,
                      double* mACsqMin, double* mACsqMax) const;
    /**
     * @brief mABsqRange_BC
     * @param mBCsq
     * @param mABsqMin
     * @param mABsqMax
     * @return
     */
    int mABsqRange_BC(const double& mBCsq,
                      double* mABsqMin, double* mABsqMax) const;

    /**
     * @brief GetLVs. Calculate LVect for mother (M) and daughter
     * particles (A, B and C) for a given Dalitz variables and under
     * the following assumptions:
     * 1. pA is directed to axis z
     * 2. xz plane is the decay plane
     * 3. pxB is ositive
     * @param mABsq. Input mAB^{2}
     * @param mACsq. Input mAC^{2}
     * @param pd. LVect of M particle (to be assigned)
     * @param pks. LVect of A particle (to be assigned)
     * @param ppip. LVect of B particle (to be assigned)
     * @param ppim. LVect of C particle (to be assigned)
     */
    void GetLVs(const double& mABsq, const double& mACsq,
                linal::LVect* pd, linal::LVect* pks,
                linal::LVect* ppip, linal::LVect* ppim) const;

    // Particle energies in the AB resonance frame
    double eA_AB(const double& mABsq) const;
    double eB_AB(const double& mABsq) const;
    double eC_AB(const double& mABsq) const;
    double eM_AB(const double& mABsq) const;

    // Particle energies in the AC resonance frame
    double eA_AC(const double& mACsq) const;
    double eB_AC(const double& mACsq) const;
    double eC_AC(const double& mACsq) const;
    double eM_AC(const double& mACsq) const;

    // Particle energies in the BC resonance frame
    double eA_BC(const double& mBCsq) const;
    double eB_BC(const double& mBCsq) const;
    double eC_BC(const double& mBCsq) const;
    double eM_BC(const double& mBCsq) const;

    // Particle momenta in the AB resonance frame
    double pA_AB(const double& mABsq) const;
    double pB_AB(const double& mABsq) const;
    double pC_AB(const double& mABsq) const;
    double pM_AB(const double& mABsq) const;

    // Particle momenta in the AC resonance frame
    double pA_AC(const double& mACsq) const;
    double pB_AC(const double& mACsq) const;
    double pC_AC(const double& mACsq) const;
    double pM_AC(const double& mACsq) const;

    // Particle momenta in the BC resonance frame
    double pA_BC(const double& mBCsq) const;
    double pB_BC(const double& mBCsq) const;
    double pC_BC(const double& mBCsq) const;
    double pM_BC(const double& mBCsq) const;

    // Resonance energy in the M frame
    double eAB_M(const double& mABsq) const;
    double eAC_M(const double& mACsq) const;
    double eBC_M(const double& mBCsq) const;

    // Resonance momentum in the M frame
    double pAB_M(const double& mABsq) const;
    double pAC_M(const double& mACsq) const;
    double pBC_M(const double& mBCsq) const;

    // * Helicity angles * //
    ///
    /// \brief CosHelAB. Cos of helisity angle for resonance R -> AB
    /// \param mABsq
    /// \param mACsq
    /// \param pq
    /// \return
    ///
    double CosHelAB(const double& mABsq, const double& mACsq,
                    double* pq) const;
    ///
    /// \brief CosHelAC. Cos of helisity angle for resonance R -> AC
    /// \param mACsq
    /// \param mBCsq
    /// \param pq
    /// \return
    ///
    double CosHelAC(const double& mACsq, const double& mBCsq,
                    double* pq) const;
    ///
    /// \brief CosHelBC. Cos of helisity angle for resonance R -> BC
    /// \param mBCsq
    /// \param mABsq
    /// \param pq
    /// \return
    ///
    double CosHelBC(const double& mBCsq, const double& mABsq,
                    double* pq) const;
    // *                                 * //

    ///
    /// \brief SetArea. Set area of the Dalitz plot
    /// \param x
    ///
    void SetArea(const double& x) {m_area = x;}
    ///
    /// \brief Area. Get area of the Dalitz plot
    /// \return
    ///
    double Area(void) const {return m_area;}

    // *** Static methods *** //
    ///
    /// \brief GetmBCsq Calculates the third Dalitz variable using
    /// other two ones
    /// \param mMsq
    /// \param mAsq
    /// \param mBsq
    /// \param mCsq
    /// \param mABsq
    /// \param mACsq
    /// \return
    ///
    static double GetmBCsq(const double& mMsq, const double& mAsq,
                           const double& mBsq, const double& mCsq,
                           const double& mABsq, const double& mACsq);
    /// Calculates momentum of either daugter in the mother particle
    /// rest frame
    static double pRes(const double& mMsq, const double& mAsq,
                       const double& mBsq);
    /// Factor for relativistic transformation between the Mo and
    /// the resonance rest frames
    static double ysq(const double& mMsq, const double& mABsq,
                      const double& mCsq);
    /// Cos of helisity angle for resonance R -> AB
    static double CosHelAB(const double& mM, const double& mA,
                           const double& mB, const double& mC,
                           const double& mABsq);
    /// Cos of helisity angle for resonance R -> AB
    static double CosHelAB(const double& mMsq, const double& mAsq,
                           const double& mBsq, const double& mCsq,
                           const double& mABsq, const double& mACsq,
                           double* pq);
    /// Cos of helisity angle for resonance R -> AC
    static double CosHelAC(const double& mMsq, const double& mAsq,
                           const double& mBsq, const double& mCsq,
                           const double& mACsq, const double& mBCsq,
                           double* pq);
    /// Cos of helisity angle for resonance R -> BC
    static double CosHelBC(const double& mMsq, const double& mAsq,
                           const double& mBsq, const double& mCsq,
                           const double& mBCsq, const double& mABsq,
                           double* pq);
    /// Min value of AB squared invariant mass
    static inline double GetmABsqMin(const double& mA, const double& mB);
    /// Max value of AB squared invariant mass
    static inline double GetmABsqMax(const double& mM, const double& mC);
    /// Momentum of daughter particle in the resonance rest frame
    static double q(const double& mRsq, const double& mA, const double& mB);
    /// Momentum of mother particle in the resonance rest frame
    static double p(const double& mRsq, const double& mM, const double& mC);
    static double eA_AB(const double& mABsq, const double& mAsq,
                        const double& mBsq);
    static double eB_AB(const double& mABsq, const double& mAsq,
                        const double& mBsq);
    static double eC_AB(const double& mABsq, const double& mMsq,
                        const double& mCsq);
    static double eM_AB(const double& mABsq, const double& mMsq,
                        const double& mCsq);
    static double eA_AC(const double& mACsq, const double& mAsq,
                        const double& mCsq);
    static double eB_AC(const double& mACsq, const double& mMsq,
                        const double& mBsq);
    static double eC_AC(const double& mACsq, const double& mAsq,
                        const double& mCsq);
    static double eM_AC(const double& mACsq, const double& mMsq,
                        const double& mBsq);
    static double eA_BC(const double& mBCsq, const double& mMsq,
                        const double& mAsq);
    static double eB_BC(const double& mBCsq, const double& mBsq,
                        const double& mCsq);
    static double eC_BC(const double& mBCsq, const double& mBsq,
                        const double& mCsq);
    static double eM_BC(const double& mBCsq, const double& mMsq,
                        const double& mAsq);

    static double pA_AB(const double& mABsq, const double& mAsq,
                        const double& mBsq);
    static double pB_AB(const double& mABsq, const double& mAsq,
                        const double& mBsq);
    static double pC_AB(const double& mABsq, const double& mMsq,
                        const double& mCsq);
    static double pM_AB(const double& mABsq, const double& mMsq,
                        const double& mCsq);
    static double pA_AC(const double& mACsq, const double& mAsq,
                        const double& mCsq);
    static double pB_AC(const double& mACsq, const double& mMsq,
                        const double& mBsq);
    static double pC_AC(const double& mACsq, const double& mAsq,
                        const double& mCsq);
    static double pM_AC(const double& mACsq, const double& mMsq,
                        const double& mBsq);

    static double pA_BC(const double& mBCsq, const double& mMsq,
                        const double& mAsq);
    static double pB_BC(const double& mBCsq, const double& mBsq,
                        const double& mCsq);
    static double pC_BC(const double& mBCsq, const double& mBsq,
                        const double& mCsq);
    static double pM_BC(const double& mBCsq, const double& mMsq,
                        const double& mAsq);

    static double eAB_M(const double& mMsq, const double& mABsq,
                        const double& mCsq);
    static double eAC_M(const double& mMsq, const double& mACsq,
                        const double& mBsq);
    static double eBC_M(const double& mMsq, const double& mBCsq,
                        const double& mAsq);

    static double pAB_M(const double& mMsq, const double& mABsq,
                        const double& mCsq);
    static double pAC_M(const double& mMsq, const double& mACsq,
                        const double& mBsq);
    static double pBC_M(const double& mMsq, const double& mBCsq,
                        const double& mAsq);

 private:
    static double ETemplateA(const double& X, const double& Y,
                             const double& Z);
    static double ETemplateC(const double& X, const double& Y,
                             const double& Z);
    static double ETemplateM(const double& X, const double& Y,
                             const double& Z);
    static double GetEnergy(const double& X, const double& Y,
                            const double& Z);
    static double Momentum(const double& E, const double& msq);
    /// Range of the (m12)^2
    static void RangeCalc(const double& e1, const double& e2,
                          const double& m1sq, const double& m2sq,
                          double* msqmin, double* msqmax);

    double m_mM;
    double m_mA;
    double m_mB;
    double m_mC;

    double m_mMsq;
    double m_mAsq;
    double m_mBsq;
    double m_mCsq;

    double m_mABsq_min, m_mABsq_max;
    double m_mACsq_min, m_mACsq_max;
    double m_mBCsq_min, m_mBCsq_max;

    double m_mass_sq_sum;
    double m_area;
};

#endif  // SRC_DALITZPHASESPACE_H_
