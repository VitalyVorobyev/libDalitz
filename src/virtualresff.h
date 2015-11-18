#ifndef VIRTUALRESFF_H
#define VIRTUALRESFF_H

#include "formfactor.h"

class VirtualResFF : public FormFactor{
public:
  VirtualResFF(const double& _r, const double& _p0);
  double operator()(const double& p) const;

private:
};

#endif // VIRTUALRESFF_H
