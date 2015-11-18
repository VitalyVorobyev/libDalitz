#ifndef B0TOD0PIPIMODEL_H
#define B0TOD0PIPIMODEL_H

#include "symdalitzmodel.h"
#include <vector>

class B0toD0pipiModel : public SymDalitzModel{
public:
  B0toD0pipiModel(void);
  B0toD0pipiModel(const double& mB, const double& mD, const double& mpi);
  EvtComplex Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
  int GetBin(const double& mp, const double& mm);

private:
  EvtComplex amp_BelleKuzmin(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
  void InitBelleModel(void);
};

#endif // B0TOD0PIPIMODEL_H
