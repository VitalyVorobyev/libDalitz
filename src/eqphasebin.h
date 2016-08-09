#ifndef EQPHASEBIN_H
#define EQPHASEBIN_H

#include "absdalitzmodel.h"

///
/// \brief The EqPhaseBin class. Calculates equal phase bin number for a given model.
///
class EqPhaseBin{
public:
  EqPhaseBin(const AbsDalitzModel* model,const int NBins);
  int      Bin(cdouble& mABsq, cdouble& mACsq);
  double delta(cdouble& mABsq, cdouble& mACsq);

private:
  double lEdge(const int i);
  double rEdge(const int i);

  const AbsDalitzModel* m_model;
  int m_nbins;

  double m_del_min;
  double m_del_max;
};

#endif // EQPHASEBIN_H
