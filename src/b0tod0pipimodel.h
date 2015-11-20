#ifndef B0TOD0PIPIMODEL_H
#define B0TOD0PIPIMODEL_H

#include "symdalitzmodel.h"
#include <vector>

class B0toD0pipiModelType{
public:
  static const int Belle = 0;
  static const int LHCb  = 1;
};

/// \brief Implementation of B0 -> D0 pi+ pi- decay amplitude
/// Belle and LHCb models are available.
/// The Belle one is from A. Kuzmin et al. (Belle Collaboration) Phys. Rev. D 76, 012006 – Published 30 July 2007
/// The LHCb one is from R. Aaij et al. (LHCb Collaboration) Phys. Rev. D 92, 032002 (2015)

class B0toD0pipiModel : public SymDalitzModel{
public:
  B0toD0pipiModel(const int type = B0toD0pipiModelType::LHCb);
  B0toD0pipiModel(const double& mB, const double& mD, const double& mpi,const int type);
  EvtComplex Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
  int GetBin(const double& mp, const double& mm);

private:
  EvtComplex amp_BelleKuzmin(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
  void InitBelleModel(void);
  int m_type;
};

#endif // B0TOD0PIPIMODEL_H
