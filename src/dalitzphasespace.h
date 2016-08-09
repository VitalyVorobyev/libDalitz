#ifndef DalitzPhaseSpace_H
#define DalitzPhaseSpace_H

#include "EvtVector4R.h"

#include <vector>
#include <complex>

typedef std::complex<double> compld;
typedef const compld         ccompld;
typedef std::vector<double>  vectd;
typedef std::vector<int>     vecti;
typedef std::vector<compld>  vectcd;
typedef const double         cdouble;
typedef const int            cint;
typedef std::string          str;
typedef const str            cstr;
typedef std::vector<str>     vectstr;

///
/// \brief The DalitzPhaseSpace class describes kinematics of a three-body decay.
///
class DalitzPhaseSpace{
public:
  ///
  /// \brief DalitzPhaseSpace
  ///
  DalitzPhaseSpace(cdouble&, cdouble&, cdouble&, cdouble&);
  ///
  /// \brief DalitzPhaseSpace
  /// \param phsp
  ///
  DalitzPhaseSpace(const DalitzPhaseSpace& phsp);
  ///
  /// \brief IsInPlot
  /// \param mABsq
  /// \param mACsq
  /// \return
  ///
  bool IsInPlot(cdouble& mABsq, cdouble& mACsq) const;
  ///
  /// \brief GetmBCsq
  /// \param mABsq
  /// \param mACsq
  /// \return
  ///
  double GetmBCsq(cdouble& mABsq, cdouble& mACsq) const;
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
  ///
  /// \brief mABsqRange_AC
  /// \param mACsq
  /// \param mABsqMin
  /// \param mABsqMax
  /// \return
  ///
  int mABsqRange_AC(cdouble& mACsq, double& mABsqMin, double& mABsqMax) const;
  ///
  /// \brief mBCsqRange_AC
  /// \param mACsq
  /// \param mBCsqMin
  /// \param mBCsqMax
  /// \return
  ///
  int mBCsqRange_AC(cdouble& mACsq, double& mBCsqMin, double& mBCsqMax) const;
  ///
  /// \brief mBCsqRange_AB
  /// \param mABsq
  /// \param mBCsqMin
  /// \param mBCsqMax
  /// \return
  ///
  int mBCsqRange_AB(cdouble& mABsq, double& mBCsqMin, double& mBCsqMax) const;
  ///
  /// \brief mACsqRange_AB
  /// \param mABsq
  /// \param mACsqMin
  /// \param mACsqMax
  /// \return
  ///
  int mACsqRange_AB(cdouble& mABsq, double& mACsqMin, double& mACsqMax) const;
  ///
  /// \brief mACsqRange_BC
  /// \param mBCsq
  /// \param mACsqMin
  /// \param mACsqMax
  /// \return
  ///
  int mACsqRange_BC(cdouble& mBCsq, double& mACsqMin, double& mACsqMax) const;
  ///
  /// \brief mABsqRange_BC
  /// \param mBCsq
  /// \param mABsqMin
  /// \param mABsqMax
  /// \return
  ///
  int mABsqRange_BC(cdouble& mBCsq, double& mABsqMin, double& mABsqMax) const;

  ///
  /// \brief GetLVs // Lorentz vectorz are calculater in M rest frame.
  /// pA is directed to axis z
  /// Decay plane is xz plane
  /// pxB is chosen to be positive
  /// \param mABsq
  /// \param mACsq
  /// \param pd
  /// \param pks
  /// \param ppip
  /// \param ppim
  ///
  void GetLVs(cdouble& mABsq,cdouble& mACsq, EvtVector4R& pd, EvtVector4R& pks,EvtVector4R& ppip,EvtVector4R& ppim) const;

  // Particle energies in the AB resonance frame
  double eA_AB(cdouble& mABsq) const;
  double eB_AB(cdouble& mABsq) const;
  double eC_AB(cdouble& mABsq) const;
  double eM_AB(cdouble& mABsq) const;
  // Particle energies in the AC resonance frame
  double eA_AC(cdouble& mACsq) const;
  double eB_AC(cdouble& mACsq) const;
  double eC_AC(cdouble& mACsq) const;
  double eM_AC(cdouble& mACsq) const;
  // Particle energies in the BC resonance frame
  double eA_BC(cdouble& mBCsq) const;
  double eB_BC(cdouble& mBCsq) const;
  double eC_BC(cdouble& mBCsq) const;
  double eM_BC(cdouble& mBCsq) const;

  // Particle momenta in the AB resonance frame
  double pA_AB(cdouble& mABsq) const;
  double pB_AB(cdouble& mABsq) const;
  double pC_AB(cdouble& mABsq) const;
  double pM_AB(cdouble& mABsq) const;
  // Particle momenta in the AC resonance frame
  double pA_AC(cdouble& mACsq) const;
  double pB_AC(cdouble& mACsq) const;
  double pC_AC(cdouble& mACsq) const;
  double pM_AC(cdouble& mACsq) const;
  // Particle momenta in the BC resonance frame
  double pA_BC(cdouble& mBCsq) const;
  double pB_BC(cdouble& mBCsq) const;
  double pC_BC(cdouble& mBCsq) const;
  double pM_BC(cdouble& mBCsq) const;

  // Resonance energy in the M frame
  double eAB_M(cdouble& mABsq) const;
  double eAC_M(cdouble& mACsq) const;
  double eBC_M(cdouble& mBCsq) const;

  // Resonance momentum in the M frame
  double pAB_M(cdouble& mABsq) const;
  double pAC_M(cdouble& mACsq) const;
  double pBC_M(cdouble& mBCsq) const;

  // * Helicity angles * //
  ///
  /// \brief CosHelAB. Cos of helisity angle for resonance R -> AB
  /// \param mABsq
  /// \param mACsq
  /// \param pq
  /// \return
  ///
  double CosHelAB(cdouble& mABsq, cdouble& mACsq, double& pq) const;
  ///
  /// \brief CosHelAC. Cos of helisity angle for resonance R -> AC
  /// \param mACsq
  /// \param mBCsq
  /// \param pq
  /// \return
  ///
  double CosHelAC(cdouble& mACsq, cdouble& mBCsq, double& pq) const;
  ///
  /// \brief CosHelBC. Cos of helisity angle for resonance R -> BC
  /// \param mBCsq
  /// \param mABsq
  /// \param pq
  /// \return
  ///
  double CosHelBC(cdouble& mBCsq, cdouble& mABsq, double& pq) const;
  // *                 * //

  ///
  /// \brief SetArea. Set area of the Dalitz plot
  /// \param x
  ///
  void SetArea(cdouble& x) {m_area = x;}
  ///
  /// \brief Area. Get area of the Dalitz plot
  /// \return
  ///
  double Area(void) const {return m_area;}

  // *** Static methods *** //
  ///
  /// \brief GetmBCsq Calculates the third Dalitz variable using other two ones
  /// \param mMsq
  /// \param mAsq
  /// \param mBsq
  /// \param mCsq
  /// \param mABsq
  /// \param mACsq
  /// \return
  ///
  static double GetmBCsq(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq,cdouble& mACsq);

  /// Calculates momentum of either daugter in the mother particle rest frame
  static double pRes(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq);

  /// Factor for relativistic transformation between the Mo and the resonance rest frames
  static double ysq(cdouble& mMsq, cdouble& mABsq, cdouble& mCsq);

  /// Cos of helisity angle for resonance R -> AB
  static double CosHelAB(cdouble& mM,cdouble& mA,cdouble& mB,cdouble& mC,cdouble& mABsq);

  /// Cos of helisity angle for resonance R -> AB
  static double CosHelAB(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mABsq, cdouble& mACsq, double& pq);

  /// Cos of helisity angle for resonance R -> AC
  static double CosHelAC(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mACsq, cdouble& mBCsq, double& pq);

  /// Cos of helisity angle for resonance R -> BC
  static double CosHelBC(cdouble& mMsq, cdouble& mAsq, cdouble& mBsq, cdouble& mCsq, cdouble& mBCsq, cdouble& mABsq, double& pq);

  /// Min value of AB squared invariant mass
  static inline double GetmABsqMin(cdouble& mA, cdouble& mB);

  /// Max value of AB squared invariant mass
  static inline double GetmABsqMax(cdouble& mM, cdouble& mC);

  /// Momentum of daughter particle in the resonance rest frame
  static double q(cdouble& mRsq, cdouble& mA, cdouble& mB);

  /// Momentum of mother particle in the resonance rest frame
  static double p(cdouble& mRsq, cdouble& mM, cdouble& mC);

  static double eA_AB(cdouble& mABsq,cdouble& mAsq,cdouble& mBsq);
  static double eB_AB(cdouble& mABsq,cdouble& mAsq,cdouble& mBsq);
  static double eC_AB(cdouble& mABsq,cdouble& mMsq,cdouble& mCsq);
  static double eM_AB(cdouble& mABsq,cdouble& mMsq,cdouble& mCsq);

  static double eA_AC(cdouble& mACsq,cdouble& mAsq,cdouble& mCsq);
  static double eB_AC(cdouble& mACsq,cdouble& mMsq,cdouble& mBsq);
  static double eC_AC(cdouble& mACsq,cdouble& mAsq,cdouble& mCsq);
  static double eM_AC(cdouble& mACsq,cdouble& mMsq,cdouble& mBsq);

  static double eA_BC(cdouble& mBCsq,cdouble& mMsq,cdouble& mAsq);
  static double eB_BC(cdouble& mBCsq,cdouble& mBsq,cdouble& mCsq);
  static double eC_BC(cdouble& mBCsq,cdouble& mBsq,cdouble& mCsq);
  static double eM_BC(cdouble& mBCsq,cdouble& mMsq,cdouble& mAsq);

  static double pA_AB(cdouble& mABsq,cdouble& mAsq,cdouble& mBsq);
  static double pB_AB(cdouble& mABsq,cdouble& mAsq,cdouble& mBsq);
  static double pC_AB(cdouble& mABsq,cdouble& mMsq,cdouble& mCsq);
  static double pM_AB(cdouble& mABsq,cdouble& mMsq,cdouble& mCsq);

  static double pA_AC(cdouble& mACsq,cdouble& mAsq,cdouble& mCsq);
  static double pB_AC(cdouble& mACsq,cdouble& mMsq,cdouble& mBsq);
  static double pC_AC(cdouble& mACsq,cdouble& mAsq,cdouble& mCsq);
  static double pM_AC(cdouble& mACsq,cdouble& mMsq,cdouble& mBsq);

  static double pA_BC(cdouble& mBCsq,cdouble& mMsq,cdouble& mAsq);
  static double pB_BC(cdouble& mBCsq,cdouble& mBsq,cdouble& mCsq);
  static double pC_BC(cdouble& mBCsq,cdouble& mBsq,cdouble& mCsq);
  static double pM_BC(cdouble& mBCsq,cdouble& mMsq,cdouble& mAsq);

  static double eAB_M(cdouble& mMsq, cdouble& mABsq, cdouble& mCsq);
  static double eAC_M(cdouble& mMsq, cdouble& mACsq, cdouble& mBsq);
  static double eBC_M(cdouble& mMsq, cdouble& mBCsq, cdouble& mAsq);

  static double pAB_M(cdouble& mMsq, cdouble& mABsq, cdouble& mCsq);
  static double pAC_M(cdouble& mMsq, cdouble& mACsq, cdouble& mBsq);
  static double pBC_M(cdouble& mMsq, cdouble& mBCsq, cdouble& mAsq);

private:
  static double ETemplateA(cdouble& X,cdouble& Y,cdouble& Z);
  static double ETemplateC(cdouble& X,cdouble& Y,cdouble& Z);
  static double ETemplateM(cdouble& X,cdouble& Y,cdouble& Z);
  static double GetEnergy( cdouble& X,cdouble& Y,cdouble& Z);
  static double Momentum(cdouble& E, cdouble& msq);
  /// Range of the (m12)^2
  static void RangeCalc(cdouble& e1,cdouble& e2,cdouble& m1sq,cdouble& m2sq,double& msqmin,double& msqmax);

  double m_mM;
  double m_mA;
  double m_mB;
  double m_mC;

  double m_mMsq;
  double m_mAsq;
  double m_mBsq;
  double m_mCsq;

  double m_mABsq_min,m_mABsq_max;
  double m_mACsq_min,m_mACsq_max;
  double m_mBCsq_min,m_mBCsq_max;

  double m_mass_sq_sum;
  double m_area;
};

typedef DalitzPhaseSpace DPhSp;
typedef const DPhSp cDPhSp;

#endif // DalitzPhaseSpace_H
