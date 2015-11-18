#ifndef GSWIDTH_H
#define GSWIDTH_H

#include "absvarwidth.h"

class GSWidth : public AbsVarWidth{
public:
  GSWidth(const double& G0, const double& m, const double& p0);

  double operator()(const double& s, const double& p) const;

private:
  double m_precalc;
};

#endif // GSWIDTH_H
