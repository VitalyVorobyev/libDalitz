#ifndef MODELINTEGRAL_H
#define MODELINTEGRAL_H

#include "absdalitzmodel.h"

///
/// \brief The ModelIntegral class. Class for calculation binned Dalitz plot parameters
/// using a gived decay amplitude model.
///
class ModelIntegral{
public:
  ModelIntegral(const AbsDalitzModel* model,const int NBins = 8, const int gridsize = 1000);

  void SetGridSize(const int gsize) {m_gsize = gsize; return;}
  void SetNBins(const int nbins)    {m_nbins = nbins; return;}
  int Calculate(cstr& label, vectd& C, vectd& S, vectd& K, vectd& Kb);

private:
  int PPbarDelta(cdouble& mABsq, cdouble& mACsq, double& P, double& Pbar, double& delta);
  const AbsDalitzModel* m_model;
  int m_nbins;
  int m_gsize;
};

#endif // MODELINTEGRAL_H
