#ifndef B0TOD0PIPIMODEL_H
#define B0TOD0PIPIMODEL_H

#include "symdalitzmodel.h"
#include "consts.h"

//#include <vector>

class B0toD0pipiModelType{
public:
  static const int Belle = 0;
  static const int LHCb  = 1;
};

///
/// \brief Implementation of B0 -> D0 pi+ pi- decay amplitude
/// Belle and LHCb models are available.
/// The Belle one is from A. Kuzmin et al. (Belle Collaboration) Phys. Rev. D 76, 012006 (2007)
/// The LHCb  one is from R. Aaij   et al. (LHCb Collaboration)  Phys. Rev. D 92, 032002 (2015)
///
class B0toD0pipiModel : public SymDalitzModel{
public:
  B0toD0pipiModel(const int type = B0toD0pipiModelType::LHCb);
  B0toD0pipiModel(const double& mB, const double& mD, const double& mpi,const int type);

private:
  void InitBelleModel(void);
  void InitLHCbModel(void);
  int m_type;
};

#endif // B0TOD0PIPIMODEL_H
