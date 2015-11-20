#ifndef KspipiModel_H
#define KspipiModel_H

#include "symdalitzmodel.h"

#include <math.h>
#include <vector>

/// \brief Implementation of D0 -> Ks0 pi+ pi- decay amplitude established in A. Poluektov et al. Phys. Rev. D 81, 112002 – Published 16 June 2010

class KspipiModel : public SymDalitzModel{
public:
  KspipiModel(void);
  KspipiModel(const double& md, const double& mks, const double& mpi);
//  EvtComplex Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
private:
//  EvtComplex amp_Belle2010(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
};

#endif // KspipiModel_H
