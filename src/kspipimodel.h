#ifndef KspipiModel_H
#define KspipiModel_H

#include "symdalitzmodel.h"

#include <math.h>
#include <vector>

class KspipiModel : public SymDalitzModel{
public:
  KspipiModel(void);
  KspipiModel(const double&, const double&, const double&);
  EvtComplex Amp(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
private:
  EvtComplex amp_Belle2010(const EvtVector4R& p4_p,const EvtVector4R& moms1,const EvtVector4R& moms2,const EvtVector4R& moms3);
};

#endif // KspipiModel_H
