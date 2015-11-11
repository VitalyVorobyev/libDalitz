#ifndef DalitzPhaseSpace_H
#define DalitzPhaseSpace_H

#include <math.h>
#include <iomanip>
#include "EvtVector4R.h"
using namespace std;

class DalitzPhaseSpace{
public:
  DalitzPhaseSpace(const double&, const double&, const double&, const double&);
  bool IsInPlot(const double& mAB,const double& mAC);
  double mBC(const double& mAC,const double& mAB);

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
