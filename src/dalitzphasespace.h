#ifndef DalitzPhaseSpace_H
#define DalitzPhaseSpace_H

#include "EvtVector4R.h"

class DalitzPhaseSpace{
public:
  DalitzPhaseSpace(const double&, const double&, const double&, const double&);
  DalitzPhaseSpace(const DalitzPhaseSpace& phsp);
  bool IsInPlot(const double& mAB,const double& mAC);
  double mBCsq(const double& mAC,const double& mAB);

  double mA() const {return m_mChA;}
  double mB() const {return m_mChB;}
  double mC() const {return m_mChC;}
  double mM() const {return m_mMo;}

  double mAB_min() const {return m_mAB_min;}
  double mAB_max() const {return m_mAB_max;}

  double mAC_min() const {return m_mAC_min;}
  double mAC_max() const {return m_mAC_max;}

  double mBC_min() const {return m_mBC_min;}
  double mBC_max() const {return m_mBC_max;}

  int mAB_range(const double&,double&,double&) const;
  void GetLVs(const double& mAB,const double& mAC, EvtVector4R& pd, EvtVector4R& pks,EvtVector4R& ppip,EvtVector4R& ppim);

  // Static methods
  static double mBCsq(const double& mMo_sq, const double& mChA_sq, const double& mChB_sq, const double& mChC_sq, const double& mABsq,const double& mACsq); /// Calculates third Dalitz variable using other two ones
  static double pRes(const double& mMo_sq, const double& mChA_sq, const double& mChB_sq);/// Calculates momentum of either daugter in the mother particle rest frame
  static double ysq(const double& mMo_sq, const double& mAB_sq, const double& mChC_sq); /// Factor for relativistic transformation between the Mo and the resonance rest frames
  static double CosHelAB(const double& mMo,const double& mA,const double& mB,const double& mC,const double& mAB); /// Cos of helisity angle for resonance R -> AB
  static inline double mABsqMin(const double& mA, const double& mB) {return pow(mA+mB,2);} /// Min value of AB squared invariant mass
  static inline double mABsqMax(const double& mM, const double& mC) {return pow(mM-mC,2);} /// Max value of AB squared invariant mass
  static double q(const double& mRsq, const double& mA, const double& mB); /// Momentum of daughter particle in the resonance rest frame
  static double p(const double& mRsq, const double& mM, const double& mC); /// Momentum of mother particle in the resonance rest frame
private:
  double m_mMo;
  double m_mChA;
  double m_mChB;
  double m_mChC;

  double m_mMo_sq;
  double m_mChA_sq;
  double m_mChB_sq;
  double m_mChC_sq;

  double m_mAB_min,m_mAB_max,m_mBC_max;
  double m_mAC_min,m_mAC_max,m_mBC_min;

  // energies in AB rest frame
  inline double eA_AC(const double& mAC) const;
  inline double eB_AC(const double& mAC) const;
};

#endif // DalitzPhaseSpace_H
