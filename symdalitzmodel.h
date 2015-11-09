#ifndef SYMDALITZMODEL_H
#define SYMDALITZMODEL_H

#include "dalitzmodel.h"

class SymDalitzModel : public DalitzModel{
public:
  SymDalitzModel(const double& mmo,const double& mcha,const double& mchb,const double& delmin,const double& delmax);

  double delta(const double& mp, const double& mm);
  void PPbarDelta(const double& mp, const double& mm, double& P, double& Pbar, double& delta);
  int GetBin(const double& mp, const double& mm);
private:
  double del_min, del_max;
  int nbins;
};

#endif // SYMDALITZMODEL_H
