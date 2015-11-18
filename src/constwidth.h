#ifndef CONSTWIDTH_H
#define CONSTWIDTH_H

#include "absvarwidth.h"

class ConstWidth : public AbsVarWidth{
public:
  ConstWidth(const double& G0);

  double operator()(const double& s = 0, const double& p = 0) const {return G0();}
};

#endif // CONSTWIDTH_H
